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

using Sundance::List;
/** 
 * Solves the coupled equations
 *
 * u_xx = v
 * v_xx = x
 * u(0) = u(1) = 0
 * v(0) = v(1) = 0
 *
 * The solution is
 * v(x) = \frac{1}{6} x (x^2 - 1)
 * u(x) = \frac{1}{120} x^5 - \frac{1}{36} x^3 + \frac{7}{360} x
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
      int nx = 128;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter rightPoint = points.subset(new RightPointTest());
      CellFilter leftPoint = points.subset(new LeftPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(5), "u");
      Expr v = new UnknownFunction(new Lagrange(3), "v");
      Expr du = new TestFunction(new Lagrange(5), "du");
      Expr dv = new TestFunction(new Lagrange(3), "dv");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(10);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, 
                          (dx*du)*(dx*u) + du*v + (dx*dv)*(dx*v) + x*dv, 
                          quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, du*u + dv*v, quad)
        + EssentialBC(rightPoint, du*u + dv*v, quad);


      /* We can now set up the linear problem! */

      LinearProblem prob(mesh, eqn, bc, List(dv,du), List(v,u), vecType);


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

      Expr x2 = x*x;
      Expr x3 = x*x2;

      Expr uExact = (1.0/120.0)*x2*x3 - 1.0/36.0 * x3 + 7.0/360.0 * x;
      Expr vExact = 1.0/6.0 * x * (x2 - 1.0);

      Expr vErr = vExact - soln[0];
      Expr uErr = uExact - soln[1];
      
      Expr vErrExpr = Integral(interior, 
                              vErr*vErr,
                              new GaussianQuadrature(6));
      
      Expr uErrExpr = Integral(interior, 
                              uErr*uErr,
                              new GaussianQuadrature(10));

      FunctionalEvaluator vErrInt(mesh, vErrExpr);
      FunctionalEvaluator uErrInt(mesh, uErrExpr);

      double uErrorSq = uErrInt.evaluate();
      std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

      double vErrorSq = vErrInt.evaluate();
      std::cerr << "v error norm = " << sqrt(vErrorSq) << std::endl << std::endl;

      double tol = 1.0e-8;
      Sundance::passFailTest(sqrt(uErrorSq+vErrorSq), tol);

    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
