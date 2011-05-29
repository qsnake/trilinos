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

#include "SundanceConstantExpr.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

ConstantExpr::ConstantExpr(const double& value)
	: SpatiallyConstantExpr(), value_(value)
{}



Set<MultipleDeriv> 
ConstantExpr::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "ConstantExpr::internalFindW found" << rtn << " for order="
    << order);

  return rtn;
}

Set<MultipleDeriv> 
ConstantExpr::internalFindV(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "ConstantExpr::internalFindV is a no-op");
  return rtn;
}


Set<MultipleDeriv> 
ConstantExpr::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "ConstantExpr::internalFindC is forwarding to findR()");
  return findR(order, context);
}


bool ConstantExpr::lessThan(const ScalarExpr* other) const
{
  const ConstantExpr* c = dynamic_cast<const ConstantExpr*>(other);
  TEST_FOR_EXCEPTION(c==0, InternalError, "cast should never fail at this point");
  return value() < c->value();
}


std::ostream& ConstantExpr::toText(std::ostream& os, bool /* paren */) const 
{
	os << value();
	return os;
}


XMLObject ConstantExpr::toXML() const 
{
	XMLObject rtn("Constant");
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}


