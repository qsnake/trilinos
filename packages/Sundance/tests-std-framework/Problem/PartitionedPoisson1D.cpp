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
 * Solves the Poisson equation in 1D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})

  int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);

    int np = MPIComm::world().getNProc();

    if (np > 1)
    {
      cout << "parallel partioned poisson test INACTIVE" << std::endl;
      return 0;
    }
    else
    {

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      //      PartitionedLineMesher::classVerbosity() = VerbExtreme;

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



      /* We can now set up the linear problem! */
      ParameterList vp = *LinearProblem::defaultVerbParams();
      LinearProblem prob(mesh, eqn, bc, v, u, vecType, vp, true); 

      LinearOperator<double> A = prob.getOperator();
      Vector<double> b = prob.getRHS()[0];
      b.scale(-1.0);

      cout << "b=" << b << std::endl;

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec-ifpack.xml"));;
#else
      ParameterXMLFileReader reader("aztec-ifpack.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cout << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Vector<double> bi = b.getBlock(0);
      Vector<double> bb = b.getBlock(1);
      LinearOperator<double> Aii = A.getBlock(0,0);
      LinearOperator<double> Aib = A.getBlock(0,1);
      LinearOperator<double> Abb = A.getBlock(1,1);
      Vector<double> xi = bi.copy();
      Vector<double> xb = bb.copy();

      solver.solve(Abb, bb, xb);
      Vector<double> c = bi - Aib * xb;
      solver.solve(Aii, c, xi);

      cout << "x_i=" << std::endl << xi << std::endl;
      cout << "x_b=" << std::endl << xb << std::endl;

      Vector<double> solnVec = prob.convertToMonolithicVector(tuple(xi), tuple(xb));
      Expr soln = prob.formSolutionExpr(tuple(solnVec));

      Expr exactSoln = x*(x-2.0);

      Expr errExpr = Integral(interior, 
        pow(soln-exactSoln, 2),
        new GaussianQuadrature(4));

      Expr derivErr = dx*(exactSoln-soln);
      Expr derivErrExpr = Integral(interior, 
        pow(dx*(soln-exactSoln), 2),
        new GaussianQuadrature(2));

      Expr fluxErrExpr = Integral(leftPoint, 
        pow(dx*(soln-exactSoln), 2),
        new GaussianQuadrature(2));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cout << "error norm = " << sqrt(errorSq) << std::endl << std::endl;



      double derivErrorSq = evaluateIntegral(mesh, derivErrExpr);
      cout << "deriv error norm = " << sqrt(derivErrorSq) << std::endl << std::endl;

      double fluxErrorSq = evaluateIntegral(mesh, fluxErrExpr);
      cout << "flux error norm = " << sqrt(fluxErrorSq) << std::endl << std::endl;

      Expr exactFluxExpr = Integral(leftPoint, 
        dx*exactSoln,
        new GaussianQuadrature(2));

      Expr numFluxExpr = Integral(leftPoint, 
        dx*soln,
        new GaussianQuadrature(2));

      cout << "exact flux = " << evaluateIntegral(mesh, exactFluxExpr) << std::endl;
      cout << "computed flux = " << evaluateIntegral(mesh, numFluxExpr) << std::endl;


      double tol = 1.0e-12;
      Sundance::passFailTest(sqrt(errorSq + derivErrorSq + fluxErrorSq), tol);
    }
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
