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

#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



DiscreteFuncElementEvaluator
::DiscreteFuncElementEvaluator(const DiscreteFuncElement* expr, 
                               const EvalContext& context)
  : SubtypeEvaluator<DiscreteFuncElement>(expr, context), 
    mi_(this->sparsity()->numDerivs()),
    miToIndexMap_(),
    stringReps_()
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "initializing discrete func evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << std::endl << *(this->sparsity()));

  static Array<string> coordNames;
  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }
  std::string funcName = expr->name();

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      /* Make sure that every derivative we're to evaluate is either
      * zero-order or a spatial derivative */
      if (this->sparsity()->deriv(i).order()==0) 
        {
          mi_[i] = MultiIndex();
        }
      else 
        {
      
          TEST_FOR_EXCEPTION(!this->sparsity()->isSpatialDeriv(i), InternalError,
                             "DiscreteFuncElementEvaluator ctor found "
                             "an entry in the sparsity superset that is not "
                             "a spatial derivative. "
                             "The bad entry is " << this->sparsity()->deriv(i) 
                             << ". The superset is " 
                             << *(this->sparsity)());

          mi_[i] = this->sparsity()->multiIndex(i);
        }
      addVectorIndex(i,i);
      TEST_FOR_EXCEPTION(miToIndexMap_.containsKey(mi_[i]), InternalError,
                         "DiscreteFuncElementEvaluator ctor detected a "
                         "duplicate multiindex");

      miToIndexMap_.put(mi_[i], i);

      if (mi_[i].order()==0)
        {
          stringReps_.append(funcName);
        }
      else
        {
          int dir = mi_[i].firstOrderDirection();
          std::string deriv = "D[" + funcName + ", " + coordNames[dir] + "]";
          stringReps_.append(deriv);
        }
    }

  
}

bool DiscreteFuncElementEvaluator::hasMultiIndex(const MultiIndex& mi) const
{
  Tabs tabs;
  bool rtn = miToIndexMap_.containsKey(mi);
  SUNDANCE_VERB_MEDIUM(tabs << "checking for mi=" << mi << " for " 
                       << expr()->toString()
                       << std::endl << tabs 
                       << " sparsity " << std::endl << *(this->sparsity()));
  
  return rtn;
}

int DiscreteFuncElementEvaluator::miIndex(const MultiIndex& mi) const
{
  return miToIndexMap_.get(mi);
}


void DiscreteFuncElementEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs;
  SUNDANCE_MSG1(mgr.verb(),
    tabs << "DiscreteFuncElementEvaluator::eval: expr=" 
    << expr()->toString());

  vectorResults.resize(mi_.size());
  for (int i=0; i<mi_.size(); i++)
    {
      vectorResults[i] = mgr.popVector();
      TEST_FOR_EXCEPTION(!vectorResults[i]->isValid(), 
                         InternalError,
                         "invalid evaluation vector allocated in "
                         "DiscreteFuncElementEvaluator::internalEval()");
      vectorResults[i]->setString(stringReps_[i]);
    }
  mgr.evalDiscreteFuncElement(expr(), mi_, vectorResults);
  mgr.stack().setVecSize(vectorResults[0]->length());
  
  if (mgr.verb() > 2)
    {
      Out::os() << tabs << "results " << std::endl;
      this->sparsity()->print(Out::os(), vectorResults,
                            constantResults);
    }
  SUNDANCE_MSG1(mgr.verb(), tabs << "DiscreteFuncEvaluator::eval() done"); 
}

