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

#ifndef SUNDANCE_UNARYEVALUATOR_H
#define SUNDANCE_UNARYEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvaluatableExpr.hpp"


namespace Sundance 
{
class EvalContext;


/**
 * 
 */
template <class ExprType> class UnaryEvaluator 
  : public SubtypeEvaluator<ExprType>
{
public:
  /** */
  UnaryEvaluator(const ExprType* expr,
    const EvalContext& context)
    : SubtypeEvaluator<ExprType>(expr, context),
      argExpr_(expr->evaluatableArg()),
      argSparsitySuperset_(argExpr_->sparsitySuperset(context)),
      argEval_(argExpr_->evaluator(context))
    {
      try
      {
        Tabs tab;

        SUNDANCE_MSG3(this->verb(), tab << "UnaryEvaluator ctor: expr = " << expr->toString());

        SUNDANCE_MSG3(this->verb(), tab << "arg sparsity superset maxOrder: " 
          << argSparsitySuperset_->maxOrder());
            
        argEval_->addClient();
            
        SUNDANCE_MSG3(this->verb(), tab << "done unary evalulator ctor");
      }
      catch(std::exception& e)
      {
        TEST_FOR_EXCEPTION(true, RuntimeError, 
          "exception detected in UnaryEvaluator: expr="
          << expr->toString() << std::endl
          << "arg=" << expr->evaluatableArg()->toString() << std::endl
          << "exception=" << e.what());
      }
    }

  /** */
  virtual ~UnaryEvaluator(){;}

  /** */
  virtual void resetNumCalls() const 
    {
      argEval_->resetNumCalls();
      Evaluator::resetNumCalls();
    }

protected:

  /** */
  const RCP<SparsitySuperset>& argSparsitySuperset() const 
    {return argSparsitySuperset_;}
      
  /** */
  const EvaluatableExpr* argExpr() const {return argExpr_;}

  /** */
  const RCP<Evaluator>& argEval() const {return argEval_;}
      

  /** */
  void evalOperand(const EvalManager& mgr,
    Array<double>& argConstantResults,
    Array<RCP<EvalVector> >& argVectorResults) const 
    {
      Tabs tabs;
      SUNDANCE_MSG1(this->verb(),  tabs << "UnaryEvaluator: evaluating operand: ");
      argEval()->eval(mgr, argConstantResults, argVectorResults);
      SUNDANCE_MSG1(this->verb(),  tabs << "UnaryEvaluator: done eval operand ");
    }
private:
  const EvaluatableExpr* argExpr_;

  RCP<SparsitySuperset> argSparsitySuperset_;

  RCP<Evaluator> argEval_;
};
}

#endif
