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

/** 
 * Solves the Poisson equation in 3D (in Q1Q1 space)
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]+1) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]+1) < 1.0e-10;})
CELL_PREDICATE(FrontPointTest, {return fabs(x[2]+1) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})
CELL_PREDICATE(BackPointTest, {return fabs(x[2]-1.0) < 1.0e-10;})

REFINE_MESH_ESTIMATE(MeshRefEst , { return 0; } , {return 1;} )
MESH_DOMAIN( MeshDomain , {return true;})

REFINE_MESH_ESTIMATE(MeshRefEst1 , { \
	
	if ( (cellPos[0] > -0.99) && (cellPos[0] < 0.5) && (cellPos[1] > -0.5) && (cellPos[1] < 0.5) &&  
	     (cellPos[2] > -0.99) && (cellPos[2] < 0.5) && (cellLevel < 1)) \
                  return 1;\
                  else return 0; } , {return 1;} )
int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();
      
      /* refinement criterion object */
      RefinementClass refCl1 = new MeshRefEst1();

      /* mesh domain, dummy class */
      MeshDomainDef meshDom = new MeshDomain();
      
      MeshType meshType = new HNMeshType3D();
      MeshSource mesher = new HNMesher3D(-1.0, -1.0, -1.0 , 2.0 , 2.0 , 2.0 , 3 , 3 , 3 , meshType , refCl1 , meshDom );
      
      Mesh mesh = mesher.getMesh();
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);
      Expr h = new CellDiameterExpr();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter faces = new DimensionalCellFilter(2);

      CellFilter left = faces.subset(new LeftPointTest());
      CellFilter right = faces.subset(new RightPointTest());
      CellFilter top = faces.subset(new TopPointTest());
      CellFilter bottom = faces.subset(new BottomPointTest());
      CellFilter front = faces.subset(new FrontPointTest());
      CellFilter back = faces.subset(new BackPointTest());

      //DOFMapBase::classVerbosity() = 5;
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      //QuadratureFamily quad4 = new GaussianQuadrature(4);

    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 0);
    watchMe.setParam("discrete function evaluation", 0);
    watchMe.setParam("integration setup", 1);
    watchMe.setParam("integral transformation", 1);
    watchMe.setParam("integration", 1);
    watchMe.setParam("fill", 5);
    watchMe.setParam("evaluation", 0);

      /* Define the weak form */
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux) + 3*vx
                          + (grad*q)*(grad*p) + 3*q,
                          quad2 );
        
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(left, vx*ux + p*q, quad2 );

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, q), 
                         List(ux, p), vecType);     

      /* */
      ParameterXMLFileReader reader("bicgstab.xml");
      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);

      FieldWriter w = new VTKWriter("Poisson_3D_Q1Q1");
      w.addMesh(mesh);

      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("p", new ExprFieldWrapper(soln[1]));
      w.write();  
      
      // test Q1 space (this is more inaccurate)
      double tol = 1.0e-2;
      Expr errExpr1 = Integral(right , (soln[0] + 6.0)*(soln[0] + 6.0)  , quad2);
      FunctionalEvaluator errInt1(mesh, errExpr1);
      double errorSq = errInt1.evaluate();
      double err_i = std::sqrt(errorSq);
      Out::os() << " error Q1 =" << err_i << std::endl;
      Sundance::passFailTest( err_i , tol);
      
      // test Q1 space (this is more inaccurate)
      tol = 1.0e-2;
      Expr errExpr = Integral( right , (soln[1] + 6.0)*(soln[1] + 6.0)  , quad2);
      FunctionalEvaluator errInt(mesh, errExpr);
      errorSq = errInt.evaluate();
      err_i = std::sqrt(errorSq);
      Out::os() << " error Q1 =" << err_i << std::endl;
      Sundance::passFailTest( err_i , tol);
   
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
