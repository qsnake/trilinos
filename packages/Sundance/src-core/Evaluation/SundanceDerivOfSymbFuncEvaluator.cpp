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

#include "SundanceDerivOfSymbFuncEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceDerivOfSymbFunc.hpp"

#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceZeroExpr.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





DerivOfSymbFuncEvaluator
::DerivOfSymbFuncEvaluator(const DerivOfSymbFunc* expr,
                           const EvalContext& context)
  : UnaryEvaluator<DerivOfSymbFunc>(expr, context),
    funcEvaluator_(),
    funcMiIndex_(-1),
    evalPtIsZero_(false),
    constResultIndex_(-1),
    funcSparsitySuperset_()
{
  Tabs tabs;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG1(verb, tabs << "initializing DerivOfSymbFunc evaluator for " 
                    << expr->toString());

  SUNDANCE_MSG2(verb, tabs << "return sparsity " << std::endl << *(this->sparsity)());

  const MultiIndex& mi = expr->mi();
  const SymbolicFuncElement* sf 
    = dynamic_cast<const SymbolicFuncElement*>(expr->evaluatableArg());

  TEST_FOR_EXCEPTION(sf==0, InternalError,
                     "Non-symbolic function detected where a symbolic "
                     "function was expected in "
                     "DerivOfSymbFuncEvaluator ctor");
  
  const TestFuncElement* t = dynamic_cast<const TestFuncElement*>(sf);
  const UnknownFuncElement* u = dynamic_cast<const UnknownFuncElement*>(sf);

  /* If we have a test function, its value will be zero. The only nonzero
   * functional derivative will be wrt D[f, mi]. The value of that deriv
   * is one. */
  if (t != 0)
    {
      constResultIndex_ = 0;
      addConstantIndex(0, constResultIndex_);
      evalPtIsZero_ = true;
      return;
    }

  /* If we have an unknown function, we need to see if it is to
   * be evaluated at zero. */
  
  TEST_FOR_EXCEPTION(u==0, InternalError,
                     "Non-unknown function detected where an unknown "
                     "function was expected in "
                     "DerivOfSymbFuncEvaluator ctor");
  
  const EvaluatableExpr* evalPt = u->evalPt();
  const ZeroExpr* z = dynamic_cast<const ZeroExpr*>(evalPt);

  if (z != 0) 
    {
      constResultIndex_ = 0;
      addConstantIndex(0, constResultIndex_);
      evalPtIsZero_ = true;
      return;
    }
  else
    {
      SUNDANCE_MSG2(verb, tabs << "setting up evaluation for discrete eval pt");
      u->setupEval(context);
    }
  
  int vecResultIndex = 0;
  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      if (this->sparsity()->state(i)==ConstantDeriv) 
        {
          constResultIndex_ = 0;
          addConstantIndex(i, constResultIndex_);
          continue;
        }
      else
        {
          addVectorIndex(i, vecResultIndex++);
          /* we have a nonzero evaluation point. We now need to find 
           * the index of my multiindex in the function's result set */

          const DiscreteFuncElement* df 
            = dynamic_cast<const DiscreteFuncElement*>(evalPt);
          
          TEST_FOR_EXCEPTION(df==0, InternalError,
                             "DerivOfSymbFuncEvaluator ctor: evaluation point of "
                             "unknown function " << u->toString() 
                             << " is not a discrete function");
  
          const SymbolicFuncElementEvaluator* uEval 
            = dynamic_cast<const SymbolicFuncElementEvaluator*>(u->evaluator(context).get());
          TEST_FOR_EXCEPTION(uEval==0, InternalError,
                             "DerivOfSymbFuncEvaluator ctor: null evaluator for "
                             "evaluation point");
          
          const DiscreteFuncElementEvaluator* dfEval = uEval->dfEval();

          TEST_FOR_EXCEPTION(dfEval==0, InternalError,
                             "DerivOfSymbFuncEvaluator ctor: evaluator for "
                             "evaluation point is not a "
                             "DiscreteFuncElementEvaluator");

          funcEvaluator_.append(dfEval);
          funcSparsitySuperset_ = dfEval->sparsity();

  
          TEST_FOR_EXCEPTION(!dfEval->hasMultiIndex(mi), InternalError,
                             "DerivOfSymbFuncEvaluator ctor: evaluator for "
                             "discrete function " << df->toString()
                             << " does not know about multiindex "
                             << mi.toString());
  
          funcMiIndex_ = dfEval->miIndex(mi);
        }
    }
}



void DerivOfSymbFuncEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults)  const
{
  Tabs tabs;
  SUNDANCE_MSG1(mgr.verb(), tabs << "DerivOfSymbFuncEvaluator::eval() expr=" 
    << expr()->toString());

  /* if the function is a test or zero func, we can take a quick way out */
  if (evalPtIsZero_)
    {
      Tabs tabs1;
      SUNDANCE_MSG2(mgr.verb(), tabs1 << "taking zero function shortcut");
      constantResults.resize(1);
      constantResults[0] = 1.0;
      return;
    }

  if (constResultIndex_ >= 0)
    {
      constantResults.resize(1);
      constantResults[0] = 1.0;
    }
  if (funcMiIndex_ >= 0)
    {
      Tabs tabs1;
      SUNDANCE_MSG2(mgr.verb(), tabs1 << "evaluating argument");
      /* evaluate the argument */
      Array<RCP<EvalVector> > funcVectorResults;
      Array<double> funcConstantResults;
      
      funcEvaluator_[0]->eval(mgr, 
        funcConstantResults, funcVectorResults);
      if (mgr.verb() > 1)
        {
          Out::os() << tabs1 << "DerivOfSymbFunc results" << std::endl;
          funcSparsitySuperset_->print(Out::os(), funcVectorResults,
                                       funcConstantResults);
        }
      vectorResults.resize(1);
      vectorResults[0] = funcVectorResults[funcMiIndex_];
    }
}



void DerivOfSymbFuncEvaluator::resetNumCalls() const 
{
  if (funcEvaluator_.size() > 0)
    {
      funcEvaluator_[0]->resetNumCalls();
    }
  Evaluator::resetNumCalls();
}
