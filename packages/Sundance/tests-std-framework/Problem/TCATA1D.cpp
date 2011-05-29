/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceParameter.hpp"
#include "SundanceExpr.hpp"


#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"

/** 
 * Solves the groundwater flow system in 1D using the NOX solver. 
 */

// column of length 1000m... hardcoded

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1000.0) < 1.0e-10;}

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      //  DOFMapBase::classVerbosity() = VerbExtreme;

      const double density = 1000.0; // kg/m^3
      const double porosity = 0.442; // dimensionless %
      const double A = 175.5; // dimensionless fit parameter
      const double B = 1.83;  // dimensionless fit parameter
      const double criticalRe = 36.73;  // dimensionless fit parameter
      const double dynvisc = 1.31;  // kg/(m-s)
      const double graindia = 1.9996e-4;  // m 
      const double charvel = 1.0;  // m/s

      
      double Reynolds = density*graindia*charvel/(dynvisc*porosity);

      Expr Re = new Parameter(Reynolds);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1000.0, 100*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellPredicate leftPointFunc = new PositionalCellPredicate(leftPointTest);
      CellPredicate rightPointFunc = new PositionalCellPredicate(rightPointTest);
      CellFilter leftPoint = points.subset(leftPointFunc);
      CellFilter rightPoint = points.subset(rightPointFunc);
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      
      Expr p = new UnknownFunction(new Lagrange(2), "p");
      Expr q = new UnknownFunction(new Lagrange(2), "q");
 
      Expr u = new TestFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      /* Define the weak form */
      Expr MassEqn = Integral(interior, q*(dx*u), quad)
	+ Integral(leftPoint, - q*u,quad)
	+ Integral(rightPoint,  - q*u,quad);
      Expr MomEqn = Integral(interior, (density/porosity)*q*q*(dx*v) + porosity*p*(dx*v) - porosity*q*v*A - (porosity*q*v*B*Re*Re)/((Re+criticalRe)*(1-porosity)), quad)
	+ Integral(leftPoint, - density*q*q*v/porosity - porosity*p*v,quad)
	+ Integral(rightPoint,- density*q*q*v/porosity - porosity*p*v,quad);

      /* Define the Dirichlet BC */
      Expr leftbc = EssentialBC(leftPoint, v*(q-charvel), quad);
      Expr rightbc = EssentialBC(rightPoint, v*(q-charvel), quad);

      /* Create a discrete space, and discretize the function 1.0 on it */
      BasisFamily L2 = new Lagrange(2);
      Array<BasisFamily> basis = tuple(L2, L2);
      DiscreteSpace discSpace(mesh, basis, vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");
      Expr p0 = u0[0];
      Expr q0 = u0[1];
     
 
/* Create a TSF NonlinearOperator object */
      std::cerr << "about to make nonlinear object" << std::endl;
      std::cerr.flush();

      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, MassEqn+MomEqn, leftbc+rightbc, Sundance::List(u,v),Sundance::List(p,q) , u0, vecType);
    
      //      F.verbosity() = VerbExtreme;
      /* Get the initial guess */
  
      Vector<double> x0 = F.getInitialGuess();
   
      
      /* Create an Aztec solver for solving the linear subproblems */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;
      
      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azOptions[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-13;
      
      LinearSolver<double> linSolver = new AztecSolver(azOptions,azParams);

      /* Now let's create a NOX solver */

      NOX::TSF::Group grp(x0, F, linSolver);

      grp.verbosity() = VerbExtreme;

      // Set up the status tests
      NOX::StatusTest::NormF statusTestA(grp, 1.0e-10);
      NOX::StatusTest::MaxIters statusTestB(20);
      NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

      // Create the list of solver parameters
      NOX::Parameter::List solverParameters;

      // Set the solver (this is the default)
      solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

      // Create the line search parameters sublist
      NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

      // Set the line search method
      lineSearchParameters.setParameter("Method","More'-Thuente");

      // Create the solver
      NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = solver.solve();

      // Print the answer
      cout << "\n" << "-- Parameter List From Solver --" << "\n";
      solver.getParameterList().print(cout);

      // Get the answer
      grp = solver.getSolutionGroup();

      // Print the answer
      cout << "\n" << "-- Final Solution From Solver --" << "\n";
      grp.print();

      

      double tol = 1.0e-12;
      Sundance::passFailTest(0, tol);

    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
