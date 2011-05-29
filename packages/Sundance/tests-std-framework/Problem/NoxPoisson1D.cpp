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
#include "NOX_TSF_StatusTestBuilder.H"

#include "TSFNOXSolver.H"
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"

/** 
 * Solves the Poisson equation in 2D using the NOX solver. 
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
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());


      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, -(dx*v)*(dx*u) - 2.0*v, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*u, quad);

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");


      /* Create a TSF NonlinearOperator object */
      NonlinearProblem prob(mesh, eqn, bc, v, u, u0, vecType);
      

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif
      ParameterList noxParams = reader.getParameters();

      std::cerr << "solver params = " << noxParams << std::endl;

      NOXSolver solver(noxParams);

      prob.solve(solver);

      Expr exactSoln = x*(x-2.0);

      Expr errExpr = Integral(interior, 
                              pow(u0-exactSoln, 2),
                              new GaussianQuadrature(4));

      Expr derivErrExpr = Integral(interior, 
                                   pow(dx*(u0-exactSoln), 2),
                                   new GaussianQuadrature(2));

      double errorSq = evaluateIntegral(mesh, errExpr);
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      double derivErrorSq = evaluateIntegral(mesh, derivErrExpr);
      std::cerr << "deriv error norm = " << sqrt(derivErrorSq) << std::endl << std::endl;

      double tol = 1.0e-12;
      Sundance::passFailTest(errorSq, tol);
      
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}

