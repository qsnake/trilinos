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
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;




ConstantEvaluator::ConstantEvaluator(const SpatiallyConstantExpr* expr, 
                                     const EvalContext& context)
  : SubtypeEvaluator<SpatiallyConstantExpr>(expr, context)
{
  Tabs tab;
  SUNDANCE_MSG1(context.setupVerbosity(), tab << "in ConstantEvaluator ctor");
  /*
   * There is only one possible nonzero derivative of this expression: the
   * zeroth-order derivative. 
   *
   * There's nothing to do in this ctor other than running some sanity checks.
   */

  TEST_FOR_EXCEPTION(this->sparsity()->numDerivs() > 1, InternalError,
                     "ConstantEvaluator ctor found a sparsity table "
                     "without more than one entry. The bad sparsity table is "
                     << *(this->sparsity)());

  if (this->sparsity()->numDerivs() > 0)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(0);

      TEST_FOR_EXCEPTION(d.order() != 0, InternalError,
                         "ConstantEvaluator ctor found a nonzero derivative "
                         "of order greater than zero. The bad sparsity "
                         "table is " << *(this->sparsity)());
      addConstantIndex(0,0);
    }
}





void ConstantEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs;
  SUNDANCE_MSG1(mgr.verb(), tabs << "ConstantEvaluator::eval() expr="
    << expr()->toString());

  if (this->sparsity()->numDerivs() > 0)
  {
    constantResults.resize(1);
    constantResults[0] = expr()->value();
    SUNDANCE_MSG2(mgr.verb(), tabs << "result=" << constantResults[0]);
  }
  else
  {
    SUNDANCE_MSG2(mgr.verb(), tabs << "no results requested");
  }
}
