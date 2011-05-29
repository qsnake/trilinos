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
#include "SundanceCellVectorExpr.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


CellVectorEvaluator::CellVectorEvaluator(const CellVectorExpr* expr, 
  const EvalContext& context)
  : SubtypeEvaluator<CellVectorExpr>(expr, context), 
    dim_(expr->dimension()),
    isTangentVector_(expr->isTangent()),
    basisVectorIndex_(expr->basisMemberIndex()),
    componentIndex_(expr->componentIndex()),
    stringRep_(expr->toString())
{

  Tabs tabs;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG1(verb, tabs << "initializing cell vector evaluator for " 
    << expr->toString());
  SUNDANCE_MSG2(verb, tabs << "return sparsity " 
    << std::endl << *(this->sparsity)());

  TEST_FOR_EXCEPTION(this->sparsity()->numDerivs() >1 , InternalError,
    "CellVectorEvaluator ctor found a sparsity table "
    "with more than one entry. The bad sparsity table is "
    << *(this->sparsity)());

  /* 
   * There is only one possible entries in the nozeros table for a
   * cell vector expression: a zeroth derivative.
   */
  
  for (int i=0; i<this->sparsity()->numDerivs(); i++)
  {
    const MultipleDeriv& d = this->sparsity()->deriv(i);

    /* for a zeroth-order derivative, evaluate the coord expr */
    TEST_FOR_EXCEPTION(d.order() != 0, RuntimeError, 
      "Derivative " << d << " is not valid for "
      "CellVectorEvaluator");
    addVectorIndex(i, 0);
  }
  
}



void CellVectorEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs;

  SUNDANCE_MSG1(mgr.verb(), tabs << "CellVectorEvaluator::eval() expr=" 
    << expr()->toString());

  SUNDANCE_MSG3(mgr.verb(), 
    tabs << "sparsity = " << std::endl << tabs << *(this->sparsity)());
  
  vectorResults.resize(1);
  vectorResults[0] = mgr.popVector();
  SUNDANCE_MSG4(mgr.verb(), tabs << "forwarding to evaluation manager");
  mgr.evalCellVectorExpr(expr(), vectorResults[0]);
  mgr.stack().setVecSize(vectorResults[0]->length());
  vectorResults[0]->setString(stringRep_);

  if (mgr.verb() > 1)
  {
    Tabs tab1;
    Out::os() << tab1 << "results " << std::endl;
    this->sparsity()->print(Out::os(), vectorResults,
      constantResults);
  }

}

