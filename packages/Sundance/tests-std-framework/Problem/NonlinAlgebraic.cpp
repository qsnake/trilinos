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

int main(int argc, char** argv)
{
  
  try
  {
    GlobalMPISession session(&argc, &argv);
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

      
    /* Create unknown and test functions, discretized using first-order
     * Lagrange interpolants */
    Expr u = new UnknownFunction(new Lagrange(1), "u");
    Expr v = new TestFunction(new Lagrange(1), "v");

    /* Create a discrete space, and discretize the function 1.0 on it */
    DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
    Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(2);

    /* Now we set up the weak form of our equation. */
    Expr eqn = Integral(interior, v*(u*u-2.0), quad);

    /* There are no boundary conditions for this problem, so the
     * BC expression is empty */
    Expr bc;

    /* We can now set up the nonlinear problem! */
    NLOp prob(mesh, eqn, bc, v, u, u0, vecType);

    /* Set up the linear solver used in solving J*delta+b = 0 */
    ParameterList solverParams;

    solverParams.set(LinearSolverBase<double>::verbosityParam(), 0);
    solverParams.set(IterativeSolver<double>::maxitersParam(), 100);
    solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-12);

    LinearSolver<double> solver = new BICGSTABSolver<double>(solverParams);

    /* Do the nonlinear solve */
      
    Vector<double> x0 = prob.getInitialGuess();
    bool converged = false;
    for (int i=0; i<20; i++)
    {
      prob.setEvalPt(x0);
      LinearOperator<double> J = prob.getJacobian();
      Vector<double> b = prob.getFunctionValue();
      Vector<double> solnVec;
      SolverState<double> state = solver.solve(J, b, solnVec);

      Out::os() << "solver state = " << std::endl << state << std::endl;

      x0 = x0 - solnVec;
      Out::os() << "step norm = " << solnVec.norm2() << std::endl;

      if (solnVec.norm2() < 1.0e-14) 
      {
        Out::os() << "Newton's method converged!" << std::endl;
        converged = true;
        break;
      }
    }
      
    if (!converged) 
    {
      Out::os() << "FAILED TO CONVERGE!" << std::endl;
    }
    else
    {
      Out::os() << "solution is " << std::endl << x0 << std::endl;
    }

    double err = 0.0;
    for (int i=0; i<x0.space().dim(); i++) 
    {
      err += pow(x0.getElement(i) - sqrt(2.0), 2.0);
    }
    err = sqrt(err)/x0.space().dim();

    Out::os() << "error norm is " << std::endl << err << std::endl;
    double tol = 1.0e-12;
    Sundance::passFailTest(err, tol);
  }
	catch(std::exception& e)
  {
    Out::os() << e.what() << std::endl;
  }
  
}
