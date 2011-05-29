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
 * Solves the Poisson equation in 2D using high-order Lagrange interpolation
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int n = 6;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, n, np,
                                                         0.0, 1.0, n, 1,
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

      
      /* Create unknown and test functions, discretized using Lagrange interpolants */
#ifdef HAVE_FIAT
      int order = 5;
      BasisFamily basis = new FIATLagrange(order);
#else
      int order = 3;
      BasisFamily basis = new Lagrange(order);
#endif
      double p = (double) order;
      Expr u = new UnknownFunction(basis, "u");
      Expr v = new TestFunction(basis, "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2*order);

      /* Define the weak form */
      Expr alpha = sqrt(2.0);
      Expr z = x + alpha*y;
      Expr exactSoln = pow(z, p);
      Expr eqn = Integral(interior, (dx*u)*(dx*v) + (dy*u)*(dy*v)  
                          + v*p*(p-1.0)*(1.0+alpha*alpha)*pow(z,p-2), quad);

      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom+top+left+right, v*(u-exactSoln), quad);


      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec.xml"));
#else
      ParameterXMLFileReader reader("aztec.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      std::cerr << "params = " << solverParams << std::endl;
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);


      Expr soln = prob.solve(solver);

      if (soln.ptr().get() != 0)
        {
      /* compute the error */
      Expr err = exactSoln - soln;
      Expr errExpr = Integral(interior, 
                              err*err,
                              quad);


      FunctionalEvaluator errInt(mesh, errExpr);

      double errorSq = errInt.evaluate();
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      Sundance::passFailTest(sqrt(errorSq), 1.0e-11);
        }
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
