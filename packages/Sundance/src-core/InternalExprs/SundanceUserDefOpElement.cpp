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

#include "SundanceUserDefOpElement.hpp"
#include "SundanceUserDefOp.hpp"
#include "SundanceUserDefOpEvaluator.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;


UserDefOpElement::UserDefOpElement(const Array<RCP<ScalarExpr> >& args,
  const RCP<Sundance::Map<EvalContext, RCP<const UserDefOpCommonEvaluator> > >& evalMap,
  const RCP<const UserDefFunctorElement>& functorElement)
  : ExprWithChildren(args), 
    commonEvaluatorsMap_(evalMap),
    functorElement_(functorElement)
{}



void UserDefOpElement::getArgDerivIndices(const Array<int>& orders,
  Map<MultiSet<int>, int>& varArgDerivs,
  Map<MultiSet<int>, int>& constArgDerivs) const
{
  int n = functorElement()->numArgs();

  int counter = 0;

  for (int o=0; o<orders.size(); o++)
  {
    int order = orders[o];
      
    if (order==0)
    {
      varArgDerivs.put(MultiSet<int>(), counter++);
    }
    else if (order==1)
    {
      for (int i=0; i<n; i++)
      {
        varArgDerivs.put(makeMultiSet<int>(i), counter++);
      }
    }
    else if (order==2)
    {
      for (int i=0; i<n; i++)
      {
        for (int j=0; j<=i; j++)
        {
          varArgDerivs.put(makeMultiSet<int>(i,j), counter++);
        }
      }
    }
    else
    {
      TEST_FOR_EXCEPTION(order > functorElement()->maxOrder() || order < 0, RuntimeError,
        "order " << order << " not supported by functor " 
        << functorElement()->masterName());
    }
  }
}

std::ostream& UserDefOpElement::toText(std::ostream& os, bool paren) const 
{
  os << functorElement()->name() << "(";
  for (int i=0; i<numChildren(); i++)
  {
    os << child(i).toString();
    if (i < numChildren()-1) os << ",";
  }
  os << ")";
  return os;
}

XMLObject UserDefOpElement::toXML() const
{
  XMLObject rtn("UserDefOpElement");
  XMLObject args("Arguments");
  for (int i=0; i<numChildren(); i++)
  {
    args.addChild(child(i).toXML());
  }
  rtn.addChild(args);
  rtn.addAttribute("op", functorElement()->name());
  return rtn;
}


Evaluator* UserDefOpElement::createEvaluator(const EvaluatableExpr* expr,
  const EvalContext& context) const
{
  const UserDefOpElement* me = dynamic_cast<const UserDefOpElement*>(expr);
  TEST_FOR_EXCEPTION(me == 0, InternalError,
    "cast failed in UserDefOpElement::createEvaluator()");
  return new UserDefOpEvaluator(me, getCommonEvaluator(context),
    context);
}



RCP<const UserDefOpCommonEvaluator> 
UserDefOpElement::getCommonEvaluator(const EvalContext& context) const 
{
  if (!commonEvaluatorsMap_->containsKey(context))
  {
    commonEvaluatorsMap_->put(context, 
      rcp(new UserDefOpCommonEvaluator(functorElement()->master(), this, context)));
  }
  return commonEvaluatorsMap_->get(context);
}
