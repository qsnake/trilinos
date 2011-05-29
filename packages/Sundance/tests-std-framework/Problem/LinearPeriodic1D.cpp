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
#include "SundancePeriodicLineMesher.hpp"
#include "SundanceUnfoldPeriodicDF.hpp"
#include "SundancePeriodicMeshType1D.hpp"

/** 
 * Solves the equation
 * 
 * \f[ u'' + 2 u' + u = \sin 2 x \f]
 *
 * with periodic BC on the interval \f$ (0, 2\pi) \f$.
 */

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      TEST_FOR_EXCEPT(np != 1);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a periodic mesh */
      int nx = 1000;
      const double pi = 4.0*atan(1.0);
      MeshType meshType = new PeriodicMeshType1D();
      MeshSource mesher = new PeriodicLineMesher(0.0, 2.0*pi, nx, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, 
        (dx*v)*(dx*u) - 2.0*v*(dx*u) - v*u + v*sin(2*x),
                          quad);
      Expr bc ; // no explicit BC needed

      /* We can now set up the linear problem! */

      LinearProblem prob(mesh, eqn, bc, v, u, vecType);


      ParameterXMLFileReader reader("amesos.xml");
      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);


      Out::os() << "solving problem " << std::endl;
      Expr soln = prob.solve(solver);

      Expr uExact = -1.0/25.0 * (4.0*cos(2.0*x) + 3.0*sin(2.0*x));

      Expr uErr = uExact - soln;
      
      Expr uErrExpr = Integral(interior, 
                              uErr*uErr,
                              new GaussianQuadrature(6));
      
      FunctionalEvaluator uErrInt(mesh, uErrExpr);

      double uErrorSq = uErrInt.evaluate();
      std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;


      /* make sure the unfolded solution is also correct */

      Out::os() << "unfolding " << std::endl;
      Expr unfoldedSoln = unfoldPeriodicDiscreteFunction(soln);
      

      Expr ufErr = uExact - unfoldedSoln;
      
      Expr ufErrExpr = Integral(interior, 
                              ufErr*ufErr,
                              new GaussianQuadrature(6));
      
      Mesh unfoldedMesh = DiscreteFunction::discFunc(unfoldedSoln)->mesh();
      FunctionalEvaluator ufErrInt(unfoldedMesh, ufErrExpr);

      double ufErrorSq = ufErrInt.evaluate();
      std::cerr << "unfolded error norm = " << sqrt(ufErrorSq) << std::endl << std::endl;

      double tol = 1.0e-3;
      Sundance::passFailTest(sqrt(uErrorSq + ufErrorSq), tol);

    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
