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
#include "SundanceUnknownParameter.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionData.hpp"


/** 
 * Solves the steady Burgers equation in 1D with forcing
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 100*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter rightPoint = points.subset(new RightPointTest());
      CellFilter leftPoint = points.subset(new LeftPointTest());
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector projector(discSpace, 1.0);
      Expr u0 = projector.project();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      /* Parameters */
      Expr p = new UnknownParameter("p");
      Expr p0 = new Sundance::Parameter(2.0);

      /* Forcing term */
      Expr f = p * (p*x*(2.0*x*x - 3.0*x + 1.0) + 2.0);

      /* Define the weak form */
      Expr eqn = Integral(interior, (dx*u)*(dx*v) + v*u*(dx*u) - v*f, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint+rightPoint, v*u, quad);

      /* Create a TSF NonlinearOperator object */
      NonlinearProblem prob(mesh, eqn, bc, v, u, u0, 
        p, p0, vecType);

      ParameterXMLFileReader reader("nox.xml");
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);
      LinearSolver<double> linSolver = solver.linSolver();

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = prob.solve(solver);
      TEST_FOR_EXCEPTION(status != NOX::StatusTest::Converged,
        runtime_error, "solve failed");


      /* compute senstivities */
      Expr sens = prob.computeSensitivities(linSolver);

      /* Write the field in ASCII format */
      FieldWriter w = new MatlabWriter("Burgers1DSoln");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(u0));
      w.addField("sens_a", new ExprFieldWrapper(sens[0]));
      w.write();


      /* check solution */
      Expr errExpr = Integral(interior, 
        pow(u0-p0*x*(1.0-x), 2),
        new GaussianQuadrature(8));
      Expr errExprA = Integral(interior, 
        pow(sens[0]-x*(1.0-x), 2),
        new GaussianQuadrature(8));

      double errorSq0 = evaluateIntegral(mesh, errExpr);
      std::cerr << "soln error norm = " << sqrt(errorSq0) << std::endl << std::endl;

      double errorSqA = evaluateIntegral(mesh, errExprA);
      std::cerr << "sens A error norm = " << sqrt(errorSqA) << std::endl << std::endl;

      double error = sqrt(errorSq0 + errorSqA);
      
      double tol = 1.0e-4;
      Sundance::passFailTest(error, tol);

    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
