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

using Sundance::List;
/** 
 * Solves the Poisson-Boltzmann equation in 2D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-2.0) < 1.0e-10;})


int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
#ifdef BLARF  

      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int n = 10;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, n, np,
                                                         0.0, 2.0, n, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      int order = 1;
      Expr u = new UnknownFunction(new Lagrange(order), "u");
      Expr v = new TestFunction(new Lagrange(order), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      double B = 0.25;
      double C = sqrt(0.5 - B*B);
      Expr exactSoln = 2.0*log(cosh(B*x + C*y));
      Expr eqn = Integral(interior, (grad*u)*(grad*v) + v*exp(-u), quad2);

      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(left+right+bottom+top, v*(u-exactSoln), quad4);

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(order), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");


      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
      

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif
      ParameterList noxParams = reader.getParameters();

      std::cerr << "solver params = " << noxParams << std::endl;

      NOXSolver solver(noxParams, F);

      solver.solve();

      Expr errExpr = Integral(interior, 
                              pow(u0-exactSoln, 2),
                              new GaussianQuadrature(4));

      Expr derivErrExpr = Integral(interior, 
                                   pow(dx*(u0-exactSoln), 2),
                                   new GaussianQuadrature(4));

      double errorSq = evaluateIntegral(mesh, errExpr);
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      double derivErrorSq = evaluateIntegral(mesh, derivErrExpr);
      std::cerr << "deriv error norm = " << sqrt(derivErrorSq) << std::endl << std::endl;

      double tol = 1.0e-12;
      Sundance::passFailTest(errorSq, tol);
#else
      Sundance::passFailTest(0.0, 1.0);
#endif
      
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 

  
}
