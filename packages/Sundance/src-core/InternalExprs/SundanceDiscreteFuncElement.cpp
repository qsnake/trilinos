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

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"

#include "SundanceDeriv.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


DiscreteFuncElement
::DiscreteFuncElement(const RCP<DiscreteFuncDataStub>& data,
  const std::string& name,
  const std::string& suffix,
  const FunctionIdentifier& fid, int myIndex)
	: EvaluatableExpr(), 
    FuncElementBase(name, suffix, fid),
    commonData_(data),
    miSet_(),
    myIndex_(myIndex)
{}


RCP<Array<Set<MultipleDeriv> > > DiscreteFuncElement
::internalDetermineR(const EvalContext& context,
                     const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "DFE::internalDetermineR() for "
                     << toString());
  Array<Set<MultipleDeriv> > RIn = RInput;
  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
    {
      const MultiIndex& mi = *i;
      int order = mi.order();
      if (order==0) RIn[0].put(MultipleDeriv());
      if (order==1) RIn[1].put(MultipleDeriv(coordDeriv(mi)));
    }

  return EvaluatableExpr::internalDetermineR(context, RIn);
}


Set<MultipleDeriv> 
DiscreteFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "DFE::internalFindW(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;

  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  if (order==0) 
    {
      if (miSet.contains(MultiIndex())) rtn.put(MultipleDeriv());
    }
  if (order==1)
    {
      for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
        {
          const MultiIndex& mi = *i;
          int diffOrder = mi.order();
          if (diffOrder==1) 
            rtn.put(MultipleDeriv(coordDeriv(mi)));
        }
    }

  return rtn;
}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "DFE::internalFindV(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;
  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  if (order==0) 
    {
      if (miSet.contains(MultiIndex())) rtn.put(MultipleDeriv());
    }
  if (order==1)
    {
      for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
        {
          const MultiIndex& mi = *i;
          int diffOrder = mi.order();
          if (diffOrder==1) 
            rtn.put(MultipleDeriv(coordDeriv(mi)));
        }
    }
  
  rtn = rtn.intersection(findR(order, context));
  return rtn;
}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_MSG5(context.setupVerbosity(), 
    tab << "DFE::internalFindC is a no-op");
  Set<MultipleDeriv> rtn;
  return rtn;
}

void DiscreteFuncElement::addMultiIndex(const MultiIndex& newMi) const
{
  miSet_.put(newMi);
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}


bool DiscreteFuncElement::lessThan(const ScalarExpr* other) const
{
  const DiscreteFuncElement* p 
    = dynamic_cast<const DiscreteFuncElement*>(other);
  TEST_FOR_EXCEPT(p==0);

  return fid() < p->fid();
}
