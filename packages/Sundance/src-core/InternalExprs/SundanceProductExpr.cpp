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


#include "SundanceProductExpr.hpp"
#include "SundanceProductEvaluator.hpp"
#include "SundanceDeriv.hpp"

#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



ProductExpr::ProductExpr(const RCP<ScalarExpr>& left,
  const RCP<ScalarExpr>& right)
	: BinaryExpr(left, right, 1)
{}


Evaluator* ProductExpr::createEvaluator(const EvaluatableExpr* expr,
  const EvalContext& context) const
{
  return new ProductEvaluator(dynamic_cast<const ProductExpr*>(expr), context);
}

bool ProductExpr::isHungryDiffOp() const
{
  return rightScalar()->isHungryDiffOp();
}



const std::string& ProductExpr::xmlTag() const 
{
	static std::string timesStr = "Times";
	static std::string divideStr = "Divide";
	if (sign() < 0) return divideStr;
	return timesStr;
}

const std::string& ProductExpr::opChar() const 
{
	static std::string timesStr = "*";
	static std::string divideStr = "/";
	if (sign() < 0) return divideStr;
	return timesStr;
}



Set<MultiSet<int> > ProductExpr::internalFindQ_W(int order, const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab0 << "ProdExpr::internalFindQ_W(" << order << ")");

  Set<MultiSet<int> > rtn;
  if (order > 2) return rtn;

  if (order==2)
  {
    rtn.put(makeMultiSet<int>(0,1));
    return rtn;
  }


  const Set<MultipleDeriv>& wLeft 
    = leftEvaluatable()->findW(0, context);
  SUNDANCE_MSG3(verb, tab0 << "wLeft=" << wLeft);

  const Set<MultipleDeriv>& wRight
    = rightEvaluatable()->findW(0, context);
  SUNDANCE_MSG3(verb, tab0 << "wRight=" << wRight);
  
  if (order==0)
  {
    if (wLeft.size() > 0)
    {
      rtn.put(makeMultiSet<int>(0));
    }
    if (wRight.size() > 0)
    {
      rtn.put(makeMultiSet<int>(1));
    }
  }
  
  if (order==1)
  {
    if (wLeft.size() > 0) rtn.put(makeMultiSet<int>(1));
    if (wRight.size() > 0) rtn.put(makeMultiSet<int>(0));
  }
  
  return rtn;
}


Set<MultiSet<int> > ProductExpr::internalFindQ_V(int order, const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab0 << "ProdExpr::internalFindQ_V(" << order << ")");

  Set<MultiSet<int> > rtn;
  if (order > 1) return rtn;

  const Set<MultipleDeriv>& vLeft 
    = leftEvaluatable()->findV(0, context);
  const Set<MultipleDeriv>& vRight
    = rightEvaluatable()->findV(0, context);

  if (order==0)
  {
    if (vLeft.size() > 0)
    {
      rtn.put(makeMultiSet<int>(0));
    }
    if (vRight.size() > 0)
    {
      rtn.put(makeMultiSet<int>(1));
    }
  }

  if (order==1)
  {
    if (vLeft.size() > 0) rtn.put(makeMultiSet<int>(1));
    if (vRight.size() > 0) rtn.put(makeMultiSet<int>(0));
  }

  
  return rtn;
}

