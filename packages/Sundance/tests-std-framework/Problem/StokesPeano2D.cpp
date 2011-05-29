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

/** 
 * Solves the Stokes equation in 2D
 */

#ifdef HAVE_SUNDANCE_PEANO
#ifdef HAVE_SUNDANCE_PEANO_NO_2D

int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy StokesPeano2D PASSED. Use Peano to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}

#else

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]+1.0) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]+1.0) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})
CELL_PREDICATE(PeggedPointTest, {return fabs(x[1]+1.0) < 1.0e-10 
                  && fabs(x[0]+1.0) < 1.0e-10 ;})

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      //int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type PeanoMeshType2D, and will
       * be built using a PeanoMesher2D. */

      MeshType meshType = new PeanoMeshType2D();
      MeshSource mesher = new PeanoMesher2D(-1.0, -1.0,  2.0 , 2.0 , 0.1 , meshType);

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
      CellFilter points = new DimensionalCellFilter(0);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      CellFilter peg = points.subset(new PeggedPointTest());


      
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
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      double beta = 0.02;
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy)
                          + (h*h*beta)*(grad*q)*(grad*p) + q*(dx*ux+dy*uy),
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = 0.5*(1.0-y*y);
      Expr bc = EssentialBC(left, vx*(ux-uInflow) + vy*uy, quad2)
        + EssentialBC(top, vx*ux + vy*uy, quad2)
        + EssentialBC(bottom, vx*ux + vy*uy, quad2)
        + EssentialBC(peg, p*q, quad2);


         

      //#define BLAHBLAH 1
#ifdef BLAHBLAH
      verbosity<Evaluator>() = 2;
      verbosity<SparsitySuperset>() = VerbExtreme;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      verbosity<Assembler>() = VerbExtreme;
#endif

      std::cerr << "Expr with children verbosity = " << verbosity<ExprWithChildren>() << std::endl;
      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, q), 
                         List(ux, uy, p), vecType);
      std::cerr << "Expr with children verbosity = " << verbosity<ExprWithChildren>() << std::endl;


#ifdef BLAHBLAH     
      cout << "row map = " << std::endl;
      prob.rowMap(0)->print(cout);
#endif
     


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

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Stokes2dPeano");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();

      

      double tol = 1.0e-3;
      Sundance::passFailTest(sqrt(errorXSq+errorYSq), tol);

    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}

#endif
#else

int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy StokesPeano2D PASSED. Use Peano to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}

#endif
