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

#include "SundanceNLOp.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceLinearSolveDriver.hpp"
#include "TSFLinearSolverDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#endif

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace std;
using namespace TSFExtended;


static Time& nlpCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("NLOp ctor"); 
  return *rtn;
}


NLOp::NLOp() 
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(),
    discreteU0_(0)
{
  TimeMonitor timer(nlpCtorTimer());
}


NLOp::NLOp(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const VectorType<double>& vecType,
  bool partitionBCs)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(u0),
    params_(),
    discreteU0_(0)
    
{
  TimeMonitor timer(nlpCtorTimer());

  Expr unkParams;
  Expr fixedParams;
  Expr fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Expr fixedFieldValues;

  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(test.flattenSpectral()), tuple(unk.flattenSpectral()), tuple(u0),
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));

  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp ctor");

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}

NLOp::NLOp(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const Expr& params, 
  const Expr& paramValues,  
  const VectorType<double>& vecType,
  bool partitionBCs)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    J_(),
    u0_(u0),
    params_(params),
    discreteU0_(0)
    
{
  TimeMonitor timer(nlpCtorTimer());

  Expr unkParams;
  Expr fixedFields;
  Expr unkParamValues;
  Expr fixedFieldValues;

  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(
            eqn, bc, tuple(test.flattenSpectral()), 
            tuple(unk.flattenSpectral()), tuple(u0), 
            unkParams, unkParamValues,
            params, paramValues,
            tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));

  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp ctor");

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}


NLOp::NLOp(const RCP<Assembler>& assembler, 
  const Expr& u0)
  : NonlinearOperatorBase<double>(),
    assembler_(assembler),
    J_(),
    u0_(u0),
    params_(),
    discreteU0_(0)
{
  TimeMonitor timer(nlpCtorTimer());
  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp ctor");

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}

TSFExtended::Vector<double> NLOp::getInitialGuess() const 
{
  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp::getInitialGuess()");
  
  Vector<double> u0 = discreteU0_->getVector();

  return u0;
}

void NLOp::setInitialGuess(const Expr& u0New) 
{
  const DiscreteFunction* in = DiscreteFunction::discFunc(u0New);
  TEST_FOR_EXCEPT(in==0);
  setEvalPt(in->getVector().copy());
  discreteU0_->setVector(currentEvalPt());
}


LinearOperator<double> NLOp::allocateJacobian() const
{
  return assembler_->allocateMatrix();
}


LinearOperator<double> NLOp::
computeJacobianAndFunction(Vector<double>& functionValue) const
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp::jacobian()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
    "null evaluation point in "
    "NLOp::jacobian()");

  discreteU0_->setVector(currentEvalPt());

  Array<Vector<double> > mv(1);
  mv[0].acceptCopyOf(functionValue);
  assembler_->assemble(J_, mv);
  functionValue.acceptCopyOf(mv[0]);

  return J_;
}

void NLOp::
computeJacobianAndFunction(LinearOperator<double>& J,
  Vector<double>& resid) const
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp::jacobian()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
    "null evaluation point in "
    "NLOp::jacobian()");

  TEST_FOR_EXCEPTION(J.ptr().get()==0, RuntimeError,
    "null Jacobian pointer in "
    "NLOp::jacobian()");

  TEST_FOR_EXCEPTION(resid.ptr().get()==0, RuntimeError,
    "null residual pointer in "
    "NLOp::jacobian()");

  discreteU0_->setVector(currentEvalPt());

  Array<Vector<double> > mv(1);
  mv[0] = resid;

  assembler_->assemble(J, mv);

  resid.acceptCopyOf(mv[0]);

  J_ = J;
}




Vector<double> NLOp::computeFunctionValue() const 
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp::computeFunctionValue()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
    "null evaluation point in "
    "NLOp::computeFunctionValue()");

  discreteU0_->setVector(currentEvalPt());

  Vector<double> rtn = createMember(range());

  Array<Vector<double> > mv(1);
  mv[0] = rtn;

  assembler_->assemble(mv);

  rtn.acceptCopyOf(mv[0]);

  return rtn;
}



void NLOp::computeFunctionValue(Vector<double>& resid) const 
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp::computeFunctionValue()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
    "null evaluation point in "
    "NLOp::computeFunctionValue()");

  TEST_FOR_EXCEPTION(resid.ptr().get()==0, RuntimeError,
    "null residual pointer in "
    "NLOp::computeFunctionValue()");

  discreteU0_->setVector(currentEvalPt());


  Array<Vector<double> > mv(1);
  mv[0] = resid;

  assembler_->assemble(mv);

  resid.acceptCopyOf(mv[0]);
}



Expr NLOp::
computeSensitivities(const LinearSolver<double>& solver) const
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
    "null discrete function pointer in "
    "NLOp::computeSensitivities()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
    "null evaluation point in "
    "NLOp::computeSensitivities()");

  TEST_FOR_EXCEPTION(J_.ptr().get()==0, RuntimeError,
    "null Jacobian pointer in "
    "NLOp::computeSensitivities()");

  TEST_FOR_EXCEPTION(params_.ptr().get()==0, RuntimeError,
    "null parameters in NLOp::computeSensitivities()");


  discreteU0_->setVector(currentEvalPt());

  Array<Vector<double> > mv(params_.size());

  assembler_->assembleSensitivities(J_, mv);

  LinearSolveDriver driver;

  Expr sens;
  int vrb = 0;
  Array<Array<string> > names(params_.size());
  for (int i=0; i<params_.size(); i++)
  {
    names[i].resize(u0_.size());
    for (int j=0; j<u0_.size(); j++) 
      names[i][j]="sens(" + u0_[j].toString() + ", " + params_[i].toString() + ")";
    mv[i].scale(-1.0);
  }

  driver.solve(solver, J_, mv, assembler_->solutionSpace(), names, vrb, sens);
  return sens;
}
