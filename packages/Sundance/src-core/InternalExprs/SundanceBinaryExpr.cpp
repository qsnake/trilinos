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

#include "SundanceBinaryExpr.hpp"
#include "SundanceExpr.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


BinaryExpr::BinaryExpr(const RCP<ScalarExpr>& left,
  const RCP<ScalarExpr>& right, int sign)
	: ExprWithChildren(tuple(left, right)), 
    sign_(sign)
{}


bool BinaryExpr::lessThan(const ScalarExpr* other) const
{
  const BinaryExpr* b = dynamic_cast<const BinaryExpr*>(other);
  TEST_FOR_EXCEPTION(b==0, InternalError, "cast should never fail at this point");
  
  if (sign_ < b->sign_) return true;
  if (sign_ > b->sign_) return false;
  return ExprWithChildren::lessThan(other);
}

std::ostream& BinaryExpr:: toText(std::ostream& os, bool paren) const 
{
	if (Expr::showAllParens() || (paren && parenthesizeSelf())) os << "(";
	leftScalar()->toText(os, parenthesizeOperands()) ;
	os	 << opChar();
	if (leftScalar()->isHungryDiffOp())
		{
			rightScalar()->toText(os, true);
		}
	else
		{
			rightScalar()->toText(os, parenthesizeOperands() || opChar() == "-");
		}
	if (Expr::showAllParens() || (paren && parenthesizeSelf())) os << ")";

	return os;
}


XMLObject BinaryExpr::toXML() const 
{
	XMLObject rtn(xmlTag());
	rtn.addChild(left().toXML());
	rtn.addChild(right().toXML());

	return rtn;
}


