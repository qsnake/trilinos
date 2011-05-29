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

#include "SundanceLinearProblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceListExpr.hpp"
#include "TSFSolverState.hpp"
#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#endif


using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;
using namespace std;


static Time& lpCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("LinearProblem ctor"); 
  return *rtn;
}


LinearProblem::LinearProblem() 
  : assembler_(),
    A_(),
    rhs_()
{
  TimeMonitor timer(lpCtorTimer());
}



LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const VectorType<double>& vecType)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(1),
    solveDriver_()
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Expr u = unk.flattenSpectral();
  Expr v = test.flattenSpectral();

  Array<Expr> zero(u.size());
  for (int i=0; i<u.size(); i++) 
  {
    Expr z = new ZeroExpr();
    zero[i] = z;
    names_[0].append(u[i].toString());
  }

  Expr u0 = new ListExpr(zero);

  Expr unkParams;
  Expr fixedParams;
  Array<Expr> fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Array<Expr> fixedFieldValues;


  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(v), tuple(u), tuple(u0),
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        fixedFields, fixedFieldValues));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));
}


LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& unkParams, 
  const Expr& unkParamVals, 
  const VectorType<double>& vecType)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(1),
    solveDriver_()
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Expr u = unk.flattenSpectral();
  Expr v = test.flattenSpectral();
  Expr alpha = unkParams.flattenSpectral();
  Expr alpha0 = unkParamVals.flattenSpectral();
  Array<Expr> zero(u.size());
  for (int i=0; i<u.size(); i++) 
  {
    Expr z = new ZeroExpr();
    zero[i] = z;
    names_[0].append(u[i].toString());
  }

  Expr u0 = new ListExpr(zero);

  Array<Expr> fixedFields;
  Expr fixedParams;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(v), tuple(u), tuple(u0),
        alpha, alpha0,
        fixedParams, fixedParams, 
        fixedFields, fixedFields));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));
}



LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const BlockArray& test, 
  const BlockArray& unk)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(unk.size()),
    solveDriver_()
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Array<Expr> v(test.size());  
  Array<Expr> u(unk.size());
  Array<Expr> u0(unk.size());

  Array<VectorType<double> > testVecType(test.size());
  Array<VectorType<double> > unkVecType(unk.size());

  for (int i=0; i<test.size(); i++)
  {
    v[i] = test[i].expr().flattenSpectral();
    testVecType[i] = test[i].vecType();
  }

  for (int i=0; i<unk.size(); i++)
  {
    u[i] = unk[i].expr().flattenSpectral();
    unkVecType[i] = unk[i].vecType();
    Array<Expr> zero(u[i].size());
    for (int j=0; j<u[i].size(); j++) 
    {
      Expr z = new ZeroExpr();
      zero[j] = z;
      names_[i].append(u[i][j].toString());
    }
    u0[i] = new ListExpr(zero);

  }

  Expr unkParams;
  Expr fixedParams;
  Array<Expr> fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Array<Expr> fixedFieldValues;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0,
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        fixedFields, fixedFieldValues));

  assembler_ = rcp(new Assembler(mesh, eqnSet, testVecType, unkVecType, partitionBCs));
}


LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const BlockArray& test, 
  const BlockArray& unk,
  const Expr& unkParams, 
  const Expr& unkParamVals)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(unk.size()),
    solveDriver_()
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Array<Expr> v(test.size());  
  Array<Expr> u(unk.size());
  Array<Expr> u0(unk.size());
  Array<VectorType<double> > testVecType(test.size());
  Array<VectorType<double> > unkVecType(unk.size());

  for (int i=0; i<test.size(); i++)
  {
    v[i] = test[i].expr().flattenSpectral();
    testVecType[i] = test[i].vecType();
  }

  for (int i=0; i<unk.size(); i++)
  {
    u[i] = unk[i].expr().flattenSpectral();
    unkVecType[i] = unk[i].vecType();
    Array<Expr> zero(u[i].size());
    for (int j=0; j<u[i].size(); j++) 
    {
      Expr z = new ZeroExpr();
      zero[j] = z;
      names_[i].append(u[i][j].toString());
    }
    u0[i] = new ListExpr(zero);
  }

  Expr fixedParams;
  Array<Expr> fixedFields;
  Expr fixedParamValues;
  Array<Expr> fixedFieldValues;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0,
        unkParams.flattenSpectral(), 
        unkParamVals.flattenSpectral(),
        fixedParams, fixedParamValues,
        fixedFields, fixedFieldValues));

  assembler_ = rcp(new Assembler(mesh, eqnSet, testVecType, unkVecType, partitionBCs));
}

LinearProblem::LinearProblem(const RCP<Assembler>& assembler)
  : assembler_(assembler),
    A_(),
    rhs_(1),
    names_()
{  
  TimeMonitor timer(lpCtorTimer());
  const RCP<EquationSet>& eqn = assembler->eqnSet();
  names_.resize(eqn->numUnkBlocks());
  for (int i=0; i<eqn->numUnkBlocks(); i++)
  {
    for (int j=0; j<eqn->numUnks(i); j++) 
    {
      names_[i].append("u(" + Teuchos::toString(i) + ", "
        + Teuchos::toString(j) + ")");
    }
  }
}

/* Return the map from cells and functions to row indices */
const RCP<DOFMapBase>& LinearProblem::rowMap(int blockRow) const 
{return assembler_->rowMap()[blockRow];}
    
/* Return the map from cells and functions to column indices */
const RCP<DOFMapBase>& LinearProblem::colMap(int blockCol) const 
{return assembler_->colMap()[blockCol];}

/* Return the discrete space in which solutions live */
const Array<RCP<DiscreteSpace> >& LinearProblem::solnSpace() const 
{return assembler_->solutionSpace();}
    
/* Return the set of row indices marked as 
 * essential boundary conditions */
const RCP<Set<int> >& LinearProblem::bcRows(int blockRow) const 
{return assembler_->bcRows()[blockRow];}

/* Return the number of block rows in the problem  */
int LinearProblem::numBlockRows() const {return assembler_->rowMap().size();}

/* Return the number of block cols in the problem  */
int LinearProblem::numBlockCols() const {return assembler_->colMap().size();}

Array<Vector<double> > LinearProblem::getRHS() const 
{
  Tabs tab;
  SUNDANCE_MSG1(assembler_->maxWatchFlagSetting("solve control"),
    tab << "LinearProblem::solve() building vector");
  assembler_->assemble(rhs_);
  return rhs_;
}


TSFExtended::LinearOperator<double> LinearProblem::getOperator() const 
{
  Tabs tab;
  SUNDANCE_MSG1(assembler_->maxWatchFlagSetting("solve control"),
    tab << "LinearProblem::solve() building matrix and vector");
  assembler_->assemble(A_, rhs_);
  return A_;
}

Expr LinearProblem::solve(const LinearSolver<double>& solver) const 
{
  Tabs tab;
  int verb = assembler_->maxWatchFlagSetting("solve control");
  Array<Vector<double> > solnVec(rhs_.size());
  
  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  for (int i=0; i<rhs_.size(); i++)
    rhs_[i].scale(-1.0);

  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() solving system");

  Expr rtn;

  /* we're not checking the status of the solve, so failures should
   * be considered fatal */
  bool save = LinearSolveDriver::solveFailureIsFatal();
  LinearSolveDriver::solveFailureIsFatal() = true;

  solveDriver_.solve(solver, A_, rhs_, solnSpace(), names_, verb, rtn);

  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() done solving system");

  /* restore original failure-handling setting */
  LinearSolveDriver::solveFailureIsFatal()=save; 

  return rtn;
}

SolverState<double> LinearProblem
::solve(const LinearSolver<double>& solver,
  Expr& soln) const 
{
  Tabs tab;
  int verb = assembler_->maxWatchFlagSetting("solve control");

  Array<Vector<double> > solnVec(rhs_.size());
  
  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  for (int i=0; i<rhs_.size(); i++)
  {
    rhs_[i].scale(-1.0);
  }

  SUNDANCE_MSG1(verb, tab << "solving LinearProblem");
  
  return solveDriver_.solve(solver, A_, rhs_, solnSpace(), names_, verb, soln);
}



Vector<double> 
LinearProblem::convertToMonolithicVector(const Array<Vector<double> >& internalBlock,
  const Array<Vector<double> >& bcBlock) const 
{return assembler_->convertToMonolithicVector(internalBlock, bcBlock);}


Expr LinearProblem::formSolutionExpr(const Array<Vector<double> >& vec) const
{
  int verb = assembler_->maxWatchFlagSetting("solve control");
  return solveDriver_.formSolutionExpr(vec, solnSpace(), names_, verb);
}



