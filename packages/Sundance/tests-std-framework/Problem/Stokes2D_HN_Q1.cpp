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
 * Solves the Stokes equation in 2D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]+1) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]+1) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})
//CELL_PREDICATE(PeggedPointTest, {return fabs(x[1]+1.0) < 1.0e-10 
//                  && fabs(x[0]+1.0) < 1.0e-10 ;})

REFINE_MESH_ESTIMATE(MeshRefEst , { return 0; } , {return 1;} )
MESH_DOMAIN( MeshDomain , {return true;})

REFINE_MESH_ESTIMATE(MeshRefEst1 , { \

	if ( (cellPos[0] > -0.5) && (cellPos[0] < 0.5) && (cellPos[1] > -0.5)
		&& (cellPos[1] < 0.5) && (cellLevel < 1)) \
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

      MeshType meshType = new HNMeshType2D();
      MeshSource mesher = new HNMesher2D(-1.0, -1.0 , 2.0 , 2.0  , 3 , 3 , meshType , refCl1 , meshDom );


      Mesh mesh = mesher.getMesh();
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr h = new CellDiameterExpr();

//       FieldWriter wMesh = new VerboseFieldWriter();
//       wMesh.addMesh(mesh);
//       wMesh.write();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      //CellFilter points = new DimensionalCellFilter(0);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      //CellFilter peg = points.subset(new PeggedPointTest());


      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr uy = new UnknownFunction(new Lagrange(1), "u_y");
      Expr vy = new TestFunction(new Lagrange(1), "v_y");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(4);
      //QuadratureFamily quad4 = new GaussianQuadrature(4);

    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 0);
    watchMe.setParam("discrete function evaluation", 0);
    watchMe.setParam("integration setup", 1);
    watchMe.setParam("integral transformation", 1);
    watchMe.setParam("integration", 6);
    watchMe.setParam("fill", 0);
    watchMe.setParam("evaluation", 0);

      /* Define the weak form */
      double beta = 0.002;
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy)
                          + (h*h*beta)*(grad*q)*(grad*p)
			  + q*(dx*ux+dy*uy)
                          ,quad2 );
		//+ Integral(left, -2*vx ,quad2 );
        
      /* Define the Dirichlet BC */
      Expr uInflow = 0.5*((1 - y*y));
      Expr bc = EssentialBC(left, vx*(ux-uInflow) , quad2 ) + //+ vy*uy
         EssentialBC(top, vx*ux + vy*uy, quad2 )
        + EssentialBC(bottom, vx*ux + vy*uy, quad2);


#ifdef BLAHBLAH
      verbosity<Evaluator>() = VerbMedium;
      verbosity<SparsitySuperset>() = VerbExtreme;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      verbosity<Assembler>() = VerbExtreme;
#endif

      std::cerr << "Expr with children verbosity = " << verbosity<ExprWithChildren>() << std::endl;
      /* We can now set up the linear problem! */
      
      LinearProblem prob(mesh, eqn, bc, List(vx, vy , q), 
                         List(ux ,uy ,p), vecType);
	 
      std::cerr << "Expr with children verbosity = " << verbosity<ExprWithChildren>() << std::endl;


#ifdef BLAHBLAH     
      cout << "row map = " << std::endl;
      prob.rowMap(0)->print(cout);
#endif
     

      ParameterXMLFileReader reader("bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      std::cerr << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);

      FieldWriter w = new VTKWriter("Stokes2d_Q1");
      w.addMesh(mesh);
    
      Expr expr_vector(List(soln[0],soln[1]));
      
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("vel", new ExprFieldWrapper(expr_vector));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();
      

      Expr exactUx = uInflow;
      Expr exactUy = 0.0;
      Expr errX = exactUx - soln[0];
      Expr errY = exactUy - soln[1];

      Expr errXExpr = Integral(interior, 
                              errX*errX,
                              new GaussianQuadrature(6));

      Expr errYExpr = Integral(interior, 
                              errY*errY,
                              new GaussianQuadrature(6));

      FunctionalEvaluator errXInt(mesh, errXExpr);
      FunctionalEvaluator errYInt(mesh, errYExpr);

      double errorXSq = errXInt.evaluate();
      double errorYSq = errYInt.evaluate();
      std::cerr << "error norm |u_x - u_x(0)| = " << sqrt(errorXSq) << std::endl << std::endl;
      std::cerr << "error norm |u_y - u_y(0)| = " << sqrt(errorYSq) << std::endl << std::endl;
      
      double tol = 1.0e-1;
      Sundance::passFailTest(sqrt(errorXSq+errorYSq), tol); 
//*/
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); 
  return Sundance::testStatus(); 
}
