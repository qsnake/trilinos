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
#include "SundanceCoordExpr.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


CoordExprEvaluator::CoordExprEvaluator(const CoordExpr* expr, 
                                       const EvalContext& context)
  : SubtypeEvaluator<CoordExpr>(expr, context), 
    doValue_(false),
    doDeriv_(false),
    stringRep_(expr->toString())
{
  int verb = context.setupVerbosity();
  Tabs tabs;
  SUNDANCE_MSG1(verb, tabs << "initializing coord expr evaluator for " 
                    << expr->toString());
  SUNDANCE_MSG2(verb, tabs << "return sparsity " << std::endl << tabs << *(this->sparsity)());

  TEST_FOR_EXCEPTION(this->sparsity()->numDerivs() > 2, InternalError,
                     "CoordExprEvaluator ctor found a sparsity table "
                     "with more than two entries. The bad sparsity table is "
                     << *(this->sparsity)());

  /* 
   * There are only two possible entries in the nozeros table for a
   * coordinate expression: a zeroth derivative, and a first-order
   * spatial derivative in the same direction as the expr's coordinate. 
   */
  
  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(i);

      /* for a zeroth-order derivative, evaluate the coord expr */
      if (d.order()==0) 
        {
          doValue_ = true;
          addVectorIndex(i, 0);
        }
      else /* for a first-order deriv, make sure it's in the proper direction,
            * then evaluate the spatial derivative. */
        {
          TEST_FOR_EXCEPTION(!this->sparsity()->isSpatialDeriv(i), InternalError,
                             "CoordExprEvaluator ctor found an entry in the "
                             "sparsity superset that is not a spatial derivative. "
                             "The bad entry is " << this->sparsity()->deriv(i) 
                             << ". The superset is " 
                             << *(this->sparsity)());

          const MultiIndex& mi = this->sparsity()->multiIndex(i);
          
          TEST_FOR_EXCEPTION(mi.order() != 1, InternalError,
                             "CoordExprEvaluator ctor found a multiindex of "
                             "order != 1. Bad multiindex is " << mi.toString());
          
          TEST_FOR_EXCEPTION(mi[expr->dir()]!=1, InternalError,
                             "CoordExprEvaluator sparsity pattern has an "
                             "element corresponding to differentiation wrt "
                             "a coordinate direction other than that of the "
                             "coord expr's direction");
          doDeriv_ = true;
          addConstantIndex(i, 0);
        }
    }
  
}



void CoordExprEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs;

  SUNDANCE_MSG1(mgr.verb(), tabs << "CoordExprEvaluator::eval() expr=" << expr()->toString());

  SUNDANCE_MSG2(mgr.verb(), tabs << "sparsity = " << std::endl 
    << *(this->sparsity)())

  if (doValue_)
    {
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << "computing value");
      vectorResults.resize(1);
      vectorResults[0] = mgr.popVector();
      mgr.evalCoordExpr(expr(), vectorResults[0]);
      mgr.stack().setVecSize(vectorResults[0]->length());
      vectorResults[0]->setString(stringRep_);
    }
  
  if (doDeriv_)
    {
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << "computing derivative");
      constantResults.resize(1);
      constantResults[0] = 1.0;
    }

  if (mgr.verb() > 2)
    {
      Tabs tab1;
      Out::os() << tab1 << "results " << std::endl;
      this->sparsity()->print(Out::os(), vectorResults,
                            constantResults);
    }

}

