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

#include "SundanceNonlinearProblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceLinearSolveDriver.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#endif


using namespace Sundance;
using namespace Teuchos;
using namespace std;
using namespace TSFExtended;


static Time& nlpCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("NonlinearProblem ctor"); 
  return *rtn;
}


NonlinearProblem::NonlinearProblem() 
  : op_()
{
  TimeMonitor timer(nlpCtorTimer());
}


NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const VectorType<double>& vecType)
  : op_(rcp(new NLOp(mesh, eqn, bc, test, unk, u0, vecType)))
{}

NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const Expr& params, 
  const Expr& paramValues,  
  const VectorType<double>& vecType)
  : op_(rcp(new NLOp(mesh, eqn, bc, test, unk, u0, 
        params, paramValues, vecType)))
{}


NonlinearProblem::NonlinearProblem(const RCP<Assembler>& assembler, 
                                   const Expr& u0)
  : op_(rcp(new NLOp(assembler, u0)))
{}



NOX::StatusTest::StatusType
NonlinearProblem::solve(const NOXSolver& solver) const
{
  RCP<NonlinearOperatorBase<double> > op = op_;
  NonlinearOperator<double> F = op;
  Vector<double> soln;
  NOX::StatusTest::StatusType rtn = solver.solve(F, soln);
  F.setEvalPt(soln);
  return rtn;
}

