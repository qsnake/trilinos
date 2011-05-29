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


#define DISABLE
#ifndef DISABLE

#include "Sundance.hpp"
#include "SundanceCubicHermite.hpp"
#include <unistd.h>
#include <sys/types.h>
/** 
 * Solves the Poisson equation in 1D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  try
    {
      Sundance::init(&argc, &argv);
      
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();
      
      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      int nx = 1;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();
      
      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter rightPoint = points.subset(new RightPointTest());
      
      /* Create unknown and test functions, discretized using cubic Hermite
       *  interpolants */
      const int p = 3;
      Expr u = new UnknownFunction(new CubicHermite(), "u" );
      Expr v = new TestFunction(new CubicHermite(), "v" );

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(6);

      WatchFlag watchMe("watch");
      watchMe.setParam("evaluation", 5);
      watchMe.setParam("integration setup", 5);
      watchMe.setParam("integration", 5);

      /* Define the weak form */
      Expr exactSoln = pow(1.0 + x, p);
      Expr source = p*(p-1.0)*pow(1.0+x, p-2.0);
      Expr eqn = Integral(interior, (dx*v)*(dx*u) + v*source, quad, watchMe);
      /* Define the Dirichlet BC */
      Expr bc;// = EssentialBC(leftPoint, v*(u-exactSoln), quad)
//        + EssentialBC(rightPoint, v*(u-exactSoln), quad);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType); 

      std::cerr << "matrix = " << std::endl << prob.getOperator() << std::endl;
      std::cerr << "rhs = " << std::endl << prob.getRHS() << std::endl;

      TEST_FOR_EXCEPT(true);
#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
#else
      ParameterXMLFileReader reader("bicgstab.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      std::cerr << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);


      Expr soln = prob.solve(solver);


      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2),
                              quad);

      double errorSq = evaluateIntegral(mesh, errExpr);
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      double tol = 1.0e-12;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}

#else

int main()
{
  return 0;
}

#endif
