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

#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceConstantEvaluator.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceParameter.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



SymbolicFuncElementEvaluator
::SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr, 
                               const EvalContext& context)
  : SubtypeEvaluator<SymbolicFuncElement>(expr, context),
    mi_(),
    spatialDerivPtrs_(),
    onePtrs_(),
    paramValuePtrs_(),
    df_(dynamic_cast<const DiscreteFuncElement*>(expr->evalPt())),
    p_(dynamic_cast<const Parameter*>(expr->evalPt())),
    stringReps_()
{
  
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing symbolic func evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << std::endl << *(this->sparsity)());

  const ZeroExpr* z 
    = dynamic_cast<const ZeroExpr*>(expr->evalPt());
  
  TEST_FOR_EXCEPTION(z==0 && df_==0, InternalError,
                     "SymbolicFuncElementEvaluator ctor detected an "
                     "evaluation point=" << expr->toString()
                     << " that is neither zero nor a discrete "
                     "function.");


  static Array<string> coordNames;
  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }
  
  int constantCounter = 0;
  int vectorCounter = 0;

  Set<MultiIndex> miSet;
  
  for (int i=0; i<this->sparsity()->numDerivs(); i++) 
    {
      if (this->sparsity()->isSpatialDeriv(i))
        {
          /* evaluate the spatial deriv applied to the evaluation point */
          TEST_FOR_EXCEPTION(z != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "spatial derivative of a zero function. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");
          TEST_FOR_EXCEPTION(p_ != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "spatial derivative of a constant parameter. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");

          mi_.append(this->sparsity()->multiIndex(i));
          miSet.put(this->sparsity()->multiIndex(i));
          addVectorIndex(i, vectorCounter);
          spatialDerivPtrs_.append(vectorCounter++);
          int dir = this->sparsity()->multiIndex(i).firstOrderDirection();
          string deriv = "D[" + df_->name() + ", " + coordNames[dir] + "]";
          stringReps_.append(deriv);
        }
      else
        {
          TEST_FOR_EXCEPTION(this->sparsity()->deriv(i).order() > 1,
                             InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "nonzero functional derivative of order greater "
                             "than one. All such derivs should have been "
                             "identified as zero by this point. The bad "
                             "derivative is " << this->sparsity()->deriv(i)
                             << ", and the bad sparsity table is "
                             << *(this->sparsity)());

          if (this->sparsity()->deriv(i).order()==0)
            {
              TEST_FOR_EXCEPTION(z != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "zero-order derivative of a zero function. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");
              /* value of zeroth functional deriv is either 
               * a discrete function or a parameter */
              if (p_ == 0)
                {
                  addVectorIndex(i, vectorCounter);
                  spatialDerivPtrs_.append(vectorCounter++);
                }
              else
                {
                  addConstantIndex(i, constantCounter);
                  paramValuePtrs_.append(constantCounter++);
                }
              mi_.append(MultiIndex());
              miSet.put(MultiIndex());
              stringReps_.append(df_->name());
            }
          else
            {
              /* value of first functional deriv is one */
              addConstantIndex(i, constantCounter);
              onePtrs_.append(constantCounter++);
            }
        }
    }

  if (p_==0 && df_ != 0)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "setting up evaluation for discrete eval pt");
      df_->setupEval(context);
      dfEval_ = dynamic_cast<const DiscreteFuncElementEvaluator*>(df_->evaluator(context).get());
    }
  else if (p_ != 0)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "setting up evaluation for parameter eval pt");
      p_->setupEval(context);
      pEval_ = dynamic_cast<const ConstantEvaluator*>(p_->evaluator(context).get());
    }
}




void SymbolicFuncElementEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs;
  
  SUNDANCE_MSG1(mgr.verb(), 
    tabs << "SymbolicFuncElementEvaluator::eval: expr=" 
    << expr()->toString());
  SUNDANCE_MSG2(mgr.verb(), tabs << "sparsity = " 
    << std::endl << *(this->sparsity)());

  constantResults.resize(onePtrs_.size() + paramValuePtrs_.size());
  vectorResults.resize(spatialDerivPtrs_.size());

  /* Evaluate discrete functions if necessary */
  if (p_==0 && df_ != 0 && mi_.size() > 0)
    {
      for (int i=0; i<mi_.size(); i++)
        {
          vectorResults[i] = mgr.popVector();
          TEST_FOR_EXCEPTION(!vectorResults[i]->isValid(), 
                             InternalError,
                             "invalid evaluation vector allocated in "
                             "SymbolicFuncElementEvaluator::internalEval()");
          vectorResults[i]->setString(stringReps_[i]);
        }
      mgr.evalDiscreteFuncElement(df_, mi_, vectorResults);
      mgr.stack().setVecSize(vectorResults[0]->length());
    }
  if (p_!=0 && mi_.size() > 0)
    {
      Array<RCP<EvalVector> > paramVectorResults;
      Array<double> paramConstResults;
      pEval_->eval(mgr, paramConstResults, paramVectorResults);
      constantResults[paramValuePtrs_[0]] = paramConstResults[0];
    }

  /* Set the known one entries to one */
  for (int i=0; i<onePtrs_.size(); i++)
    {
      constantResults[onePtrs_[i]] = 1.0;
    }
  if (verb() > 2)
    {
      Out::os() << tabs << "results " << std::endl;
      this->sparsity()->print(Out::os(), vectorResults,
                            constantResults);
    }
  SUNDANCE_MSG1(mgr.verb(), tabs << "SymbolicFuncEvaluator::eval() done"); 

}


