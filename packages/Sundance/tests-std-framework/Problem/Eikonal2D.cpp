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

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"


#ifdef Trilinos_DATA_DIR
/** 
 * Solves the Eikonal equation in 2D
 * | \nabla u |^2 = 1 + \epsilon \nabla^2 u
 */


CELL_PREDICATE(BdryPointTest,
               {
                 double r2 = x[0]*x[0]+x[1]*x[1];
                 bool rtn = ::fabs(r2 - 100.0) < 1.0e-6;
                 return rtn;
               })

extern "C" {double besi0_(double* x);}

class BesselI0Func : public PointwiseUserDefFunctor0
{
public: 
  /** */
  BesselI0Func() : PointwiseUserDefFunctor0("I0", 1, 1) {;}

  /** */
  void eval0(const double* vars, double* f) const 
  {
    double* x = const_cast<double*>(vars);
    f[0] = besi0_(x);
  }
};

Expr BesselI0(const Expr& x)
{
  Expr rtn = new UserDefOp(x, rcp(new BesselI0Func()));
  return rtn;
}

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      double epsilon = 2.0;   

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh object to fill in from a file. It will be of type BasisSimplicialMesh. */
      MeshType meshType = new BasicSimplicialMeshType();

      /* read the mesh from the file disk.1*/
      MeshSource meshReader = new TriangleMeshReader("disk.1", meshType);

      Mesh mesh = meshReader.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain.  also define edges and boundaries.*/
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter bdry = edges.subset(new BdryPointTest);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create a discrete space, and discretize the function 1.0 on it to represent rhs */
      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      Expr S = new DiscreteFunction(discSpace, 1.0, "S");

      /* Discretize the function 0.0 on the same discrete space as an initial guess*/
      Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");

      /* Create differential operators, gradient, and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = Sundance::List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);
     
      /* Define the weak form */
      Expr eqn = Integral(interior, (dx*u)*(dx*u)*v+(dy*u)*(dy*u)*v+epsilon*(grad*v)*(grad*u)-S*v, quad2);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bdry, v*u, quad4);

      /* Create a TSF NonlinearOperator object */
      NonlinearProblem nlp(mesh, eqn, bc, v, u, u0, vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif
      ParameterList noxParams = reader.getParameters();

      std::cerr << "solver params = " << noxParams << std::endl;

      NOXSolver solver(noxParams);

      /* Solve the problem. The solution is written into the expr u0 */
      nlp.solve(solver);

     

    


      /* Check against exact solution */
      double R0 = 10.0;
      double z = R0/epsilon;
      double pi = 4.0*atan(1.0);
      double I0 = besi0_(&z);
      Expr r = sqrt(x*x + y*y);

//      Expr I0r = new UserDefOp(r/epsilon, rcp(new BesselI0Func()));
      Expr exactSoln = epsilon * (log(I0) - log(BesselI0(r/epsilon)));
//      Expr exactSoln = epsilon * (log(I0) - log(I0r));
      Expr exactDisc = L2Projector(discSpace, exactSoln).project();

      
      /* this code writes the result to a file so we can visualize using paraview */
      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Eikonal2D");
      w.addMesh(mesh);
      w.addField("numerical soln", new ExprFieldWrapper(u0[0]));
      w.addField("exact soln", new ExprFieldWrapper(exactDisc));
      w.write();

      Expr errExpr = Integral(interior, 
                              pow(u0[0]-exactSoln, 2.0)/(epsilon+pow(exactSoln,2.0)),
                              new GaussianQuadrature(4) );
      double errorSq = evaluateIntegral(mesh, errExpr)/pi/R0/R0;
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(::sqrt(errorSq), tol);
      
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}


#else



int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy Eikonal2D PASSED. Enable Trilinos_DATA_DIR to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif
