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

#include "SundanceSumExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"



using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;


SumExpr::SumExpr(const RCP<ScalarExpr>& left,
  const RCP<ScalarExpr>& right, int sign)
	: BinaryExpr(left, right, sign), sumTree_()
{
  /*
    Expr L = Expr::handle(left);
    Expr R = Expr::handle(right);

    sumTree_ = L.getSumTree();
    Map<Expr, int> rightTree = R.getSumTree();

    for (Map<Expr, int>::const_iterator i=rightTree.begin(); i!=rightTree.end(); i++)
    {
    int leftCount = 0;
    if (sumTree_.containsKey(i->first))
    {
    leftCount = sumTree_[i->first];
    }
    int rightCount = sign * i->second;
    sumTree_.put(i->first, leftCount + rightCount);
    }
  */
}

bool SumExpr::isHungryDiffOp() const
{
  return leftScalar()->isHungryDiffOp() || rightScalar()->isHungryDiffOp();
}


const std::string& SumExpr::xmlTag() const 
{
	static std::string plusStr = "Plus";
	static std::string minusStr = "Minus";
	if (sign() < 0) return minusStr;
	return plusStr;
}

const std::string& SumExpr::opChar() const 
{
	static std::string plusStr = "+";
	static std::string minusStr = "-";
	if (sign() < 0) return minusStr;
	return plusStr;
}


bool SumExpr::everyTermHasTestFunctions() const
{
  return leftEvaluatable()->everyTermHasTestFunctions()
    && rightEvaluatable()->everyTermHasTestFunctions();
}

bool SumExpr::isLinearInTests() const
{
  bool leftHasTests = leftScalar()->hasTestFunctions();
  bool rightHasTests = rightScalar()->hasTestFunctions();

  bool leftIsLinear = leftScalar()->isLinearInTests();
  bool rightIsLinear = rightScalar()->isLinearInTests();

  return (!leftHasTests || leftIsLinear) && (!rightHasTests || rightIsLinear);
}


bool SumExpr::isLinearForm(const Expr& u) const 
{
  bool LL = leftScalar()->isLinearForm(u);
  bool RL = rightScalar()->isLinearForm(u);
  bool LI = leftScalar()->isIndependentOf(u);
  bool RI = rightScalar()->isIndependentOf(u);

  return ( (LL && (RL || RI)) || (RL && (LL || LI)) );
}

bool SumExpr::isQuadraticForm(const Expr& u) const
{
  bool LQ = leftScalar()->isQuadraticForm(u);
  bool RQ = rightScalar()->isQuadraticForm(u);
  bool LL = leftScalar()->isLinearForm(u);
  bool RL = rightScalar()->isLinearForm(u);
  bool LI = leftScalar()->isIndependentOf(u);
  bool RI = rightScalar()->isIndependentOf(u);

  return ( (LQ && (RQ || RL || RI)) || (RQ && (LQ || LL || LI))); 
}
