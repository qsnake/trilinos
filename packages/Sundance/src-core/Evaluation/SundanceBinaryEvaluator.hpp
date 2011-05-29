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

#ifndef SUNDANCE_BINARYEVALUATOR_H
#define SUNDANCE_BINARYEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvalManager.hpp"

namespace Sundance 
{
class EvalContext;
  


/**
 * 
 */
template <class ExprType> class BinaryEvaluator 
  : public SubtypeEvaluator<ExprType>
{
public:
  /** */
  BinaryEvaluator(const ExprType* expr,
    const EvalContext& context)
    : SubtypeEvaluator<ExprType>(expr, context),
      leftExpr_(expr->leftEvaluatable()),
      rightExpr_(expr->rightEvaluatable()),
      leftSparsity_(leftExpr_->sparsitySuperset(context)),
      rightSparsity_(rightExpr_->sparsitySuperset(context)),
      leftEval_(leftExpr_->evaluator(context)),
      rightEval_(rightExpr_->evaluator(context))
    {
      Tabs tab;
      leftEval_->addClient();
      rightEval_->addClient();
    }

  /** */
  virtual ~BinaryEvaluator(){;}

  /** */
  virtual void resetNumCalls() const 
    {
      leftEval_->resetNumCalls();
      rightEval_->resetNumCalls();
      Evaluator::resetNumCalls();
    }

protected:
      
  /** */
  const RCP<SparsitySuperset>& leftSparsity() const 
    {return leftSparsity_;}
      
  /** */
  const RCP<SparsitySuperset>& rightSparsity() const 
    {return rightSparsity_;}

  /** */
  const EvaluatableExpr* leftExpr() const {return leftExpr_;}

  /** */
  const EvaluatableExpr* rightExpr() const {return rightExpr_;}

  /** */
  const RCP<Evaluator>& leftEval() const {return leftEval_;}

  /** */
  const RCP<Evaluator>& rightEval() const {return rightEval_;}

  /** */
  void evalChildren(const EvalManager& mgr,
    Array<double>& leftConstResults,
    Array<RCP<EvalVector> >& leftVecResults,
    Array<double>& rightConstResults,
    Array<RCP<EvalVector> >& rightVecResults) const 
    {
      Tabs tabs;
      SUNDANCE_MSG2(mgr.verb(), 
        tabs << "Evaluating left and right children: "
        << std::endl << tabs << "left=" << leftExpr_->toString()
        << std::endl << tabs << "right=" << rightExpr_->toString());
      SUNDANCE_MSG2(mgr.verb(),
        tabs << "Evaluating left=" << leftExpr_->toString());
      leftEval()->eval(mgr, leftConstResults, leftVecResults);
        
      SUNDANCE_MSG2(mgr.verb(),
        tabs << "Evaluating right=" << rightExpr_->toString());
      rightEval()->eval(mgr, rightConstResults, rightVecResults);
    }

private:
  const EvaluatableExpr* leftExpr_;

  const EvaluatableExpr* rightExpr_;

  RCP<SparsitySuperset> leftSparsity_;

  RCP<SparsitySuperset> rightSparsity_;

  RCP<Evaluator> leftEval_;

  RCP<Evaluator> rightEval_;
};
}


#endif
