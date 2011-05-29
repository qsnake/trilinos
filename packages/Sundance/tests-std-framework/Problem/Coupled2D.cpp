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
#include <unistd.h>
#include <sys/unistd.h>

using Sundance::List;
/** 
 * Solves the coupled equations
 *
 * u_xx = v
 * v_xx = 1
 * u(0) = u(1) = 0
 * v(0) = v(1) = 0
 *
 * The solution is
 * v(x) = -1/2 x (1-x)
 * u(x) = 1/24 x (x^3 - 2 x^2 + 1)
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
      int np = MPIComm::world().getNProc();
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEST_FOR_EXCEPT(npx < 1);
      TEST_FOR_EXCEPT(npy < 1);
      TEST_FOR_EXCEPT(npx * npy != np);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 32;

      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx,
                                                         0.0, 2.0, nx, npy,
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
      Expr u = new UnknownFunction(new Lagrange(3), "u");
      Expr v = new UnknownFunction(new Lagrange(2), "v");
      Expr du = new TestFunction(new Lagrange(3), "du");
      Expr dv = new TestFunction(new Lagrange(2), "dv");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);
      Expr dy = new Derivative(1);
      Expr y = new CoordExpr(1);
      Expr grad = List(dx, dy);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(6);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, 
                          (grad*du)*(grad*u) + du*v + (grad*dv)*(grad*v) + x*dv, 
                          quad2);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(left, du*u + dv*v, quad2)
        + EssentialBC(right, du*u + dv*v, quad2);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(dv,du), List(v,u), vecType);

      

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      //azOptions[AZ_ml] = 1;
      //azOptions[AZ_ml_levels] = 4;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-8;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);



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
                              new GaussianQuadrature(6));

      FunctionalEvaluator vErrInt(mesh, vErrExpr);
      FunctionalEvaluator uErrInt(mesh, uErrExpr);

      double uErrorSq = uErrInt.evaluate();
      std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

      double vErrorSq = vErrInt.evaluate();
      std::cerr << "v error norm = " << sqrt(vErrorSq) << std::endl << std::endl;

      
      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Coupled2D");
      w.addMesh(mesh);
      w.addField("v", new ExprFieldWrapper(soln[0]));
      w.addField("u", new ExprFieldWrapper(soln[1]));
      w.write();



      double tol = 1.0e-5;
      Sundance::passFailTest(sqrt(uErrorSq+vErrorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
