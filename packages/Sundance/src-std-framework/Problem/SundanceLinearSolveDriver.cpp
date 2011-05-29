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

#include "SundanceLinearSolveDriver.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearSolverImpl.hpp"
#endif




using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;
using namespace std;


SolverState<double> 
LinearSolveDriver::solve(const LinearSolver<double>& solver,
  const LinearOperator<double>& A,
  const Array<Vector<double> >& rhs,
  const Array<RCP<DiscreteSpace> >& solutionSpace,
  const Array<Array<string> >& names,
  int verb,
  Expr& soln) const
{
  Tabs tab(0);
  Array<Vector<double> > solnVec(rhs.size());
  SolverState<double> state;

  for (int i=0; i<rhs.size(); i++)
  {
    Tabs tab1;

    solnVec[i] = rhs[i].copy();
    
    SUNDANCE_MSG2(verb, tab1 << "solving with RHS #" << i 
      << " of " << rhs.size());

    state = solver.solve(A, rhs[i], solnVec[i]);
    
    SUNDANCE_MSG2(verb, tab1 << "solve completed with status="
      << state.stateDescription());

    /* deal with a failure to converge */
    if (state.finalState() != SolveConverged)
    {
      TeuchosOStringStream ss;
      ss << "Solve failed! state = "
         << state.stateDescription()
         << "\nmessage=" << state.finalMsg()
         << "\niters taken = " << state.finalIters()
         << "\nfinal residual = " << state.finalResid();

      /* If requested, write the bad matrix and vector */
      if (dumpBadMatrix())
      {
        if (A.ptr().get() != 0)
        {
          ofstream osA(badMatrixFilename().c_str());
          A.print(osA);
          ss << "\nmatrix written to " << badMatrixFilename();
        }
        else
        {
          ss << "\nthe matrix is null! Evil is afoot in your code...";
        }
        if (rhs[i].ptr().get() != 0)
        {
          ofstream osb(badVectorFilename().c_str());
          rhs[i].print(osb);
          ss << "\nRHS vector written to " << badVectorFilename();
        }
        else
        {
          ss << "\nthe RHS vector is null! Evil is afoot in your code...";
        }
      }
      
      /* If solve errors are fatal, throw an exception */
      TEST_FOR_EXCEPTION(solveFailureIsFatal(),
        RuntimeError, TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss));

      /* otherwise, return the state information */
      return state;
    }
  }
   
  /* Put the solution vector into a discrete function */
  if (soln.ptr().get()==0)
  {
    soln = formSolutionExpr(solnVec, solutionSpace, names, verb);
  }
  else
  {
    writeIntoSolutionExpr(solnVec, soln, verb);
  }

  return state;
      
}



Expr LinearSolveDriver::formSolutionExpr(
  const Array<Vector<double> >& solnVector,
  const Array<RCP<DiscreteSpace> >& solutionSpace,
  const Array<Array<string> >& names,
  int verb) const
{
  Array<Expr> cols(solnVector.size());

  for (int m=0; m<solnVector.size(); m++)
  {
    Array<Expr> col(solutionSpace.size());

    for (int i=0; i<col.size(); i++)
    {
      std::string name = names[m][i];
      if (col.size() > 1) name += "[" + Teuchos::toString(i) + "]";
      col[i] = new DiscreteFunction(*(solutionSpace[i]),
        solnVector[m].getBlock(i), name);
    }
    if (col.size() > 1)
    {
      cols[m] = new ListExpr(col);
    }
    else
    {
      cols[m] = col[0];
    }
  }

  if (cols.size() > 1)
  {
    return new ListExpr(cols);;
  }
  else
  {
    return cols[0];
  }
}


void LinearSolveDriver::writeIntoSolutionExpr(
  const Array<Vector<double> >& solnVector,
  Expr soln,
  int verb) const 
{
  TEST_FOR_EXCEPT(soln.size() != solnVector.size());
  for (int i=0; i<solnVector.size(); i++)
  {
    Expr u = soln[i];
    DiscreteFunction::discFunc(u)->setVector(solnVector[i]);
  }
}
