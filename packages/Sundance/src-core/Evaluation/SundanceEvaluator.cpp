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

#include "SundanceEvaluator.hpp"
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




Evaluator::Evaluator()
  : numClients_(0),
    numCalls_(0),
    vectorResultCache_(),
    constantResultCache_(),
    constantIndexMap_(),
    vectorIndexMap_(),
    vectorIndices_(),
    constantIndices_()
{}



void Evaluator::eval(const EvalManager& mgr,
                     Array<double>& constantResults,
                     Array<RCP<EvalVector> >& vectorResults) const
{
  Tabs tabs;

  if (numCalls_ == 0)
    {
      internalEval(mgr, constantResultCache_, vectorResultCache_);
    }
  
  numCalls_++;

  /* Go ahead and copy the constant results every time, 
   * since this is cheap */
  if (constantResultCache_.size() > 0) constantResults = constantResultCache_;

  /* If all clients have called, we can return the original data
   * which can then be changed by the client. */
  if (numCalls_ == numClients_)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "surrendering cached results");
      vectorResults = vectorResultCache_;
    }
  else /* Otherwise, make a copy of the data */
    {
      SUNDANCE_VERB_MEDIUM(tabs << "cloning cached results");
      vectorResults.resize(vectorResultCache_.size());
      for (int i=0; i < vectorResults.size(); i++)
        {
          vectorResults[i] = vectorResultCache_[i]->clone();
        }
    }
}


void Evaluator::addConstantIndex(int index, int constantIndex)
{
  TEST_FOR_EXCEPTION(constantIndexMap_.containsKey(index), InternalError,
                     "duplicate index " << index 
                     << " found in Evaluator::addConstantIndex");
  constantIndexMap_.put(index, constantIndex);
  constantIndices_.append(index);
}

void Evaluator::addVectorIndex(int index, int vectorIndex)
{
  TEST_FOR_EXCEPTION(vectorIndexMap_.containsKey(index), InternalError,
                     "duplicate index " << index 
                     << " found in Evaluator::addVectorIndex");
  vectorIndexMap_.put(index, vectorIndex);
  vectorIndices_.append(index);
}

