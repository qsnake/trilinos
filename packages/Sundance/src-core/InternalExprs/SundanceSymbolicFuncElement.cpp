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

#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceZeroExpr.hpp"

#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



SymbolicFuncElement::SymbolicFuncElement(const std::string& name,
  const std::string& suffix,
  const FunctionIdentifier& fid,
  const RCP<const CommonFuncDataStub>& data)
	: EvaluatableExpr(), FuncElementBase(name, suffix, fid),
    commonData_(data),
    evalPt_(),
    evalPtDerivSetIndices_()
{}

void SymbolicFuncElement::registerSpatialDerivs(const EvalContext& context, 
                                                const Set<MultiIndex>& miSet) const
{
  evalPt()->registerSpatialDerivs(context, miSet);
  EvaluatableExpr::registerSpatialDerivs(context, miSet);
}

void SymbolicFuncElement::accumulateFuncSet(Set<int>& funcDofIDs, 
  const Set<int>& activeSet) const
{
  if (activeSet.contains(fid().dofID())) 
  {
    funcDofIDs.put(fid().dofID());
  }
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "SFE::internalFindW(order=" << order << ") for "
                     << toString());

  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "findW() for eval point");
    evalPt()->findW(order, context);
  }

  if (order==0) 
    {
      Tabs tab1;
      if (!evalPtIsZero()) 
        {
          SUNDANCE_MSG5(verb, tab1 << "value of " << toString() << " is nonzero" );
          rtn.put(MultipleDeriv());
        }
      else
        {
          SUNDANCE_MSG5(verb,  tab1 << "value of " << toString() << " is zero" );
        }
    }
  else if (order==1)
    {
      Deriv d = funcDeriv(this);
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }
  
  SUNDANCE_MSG3(verb,  tab << "SFE: W[" << order << "] = " << rtn );

  return rtn;
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "SFE::internalFindV(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "findV() for eval point");
    evalPt()->findV(order, context);
  }

  if (order==0) 
    {
      if (!evalPtIsZero()) rtn.put(MultipleDeriv());
    }

  SUNDANCE_MSG5(verb,  tab << "SFE: V = " << rtn );
  SUNDANCE_MSG5(verb,  tab << "SFE: R = " << findR(order, context) );
  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_MSG3(verb,  tab << "SFE: V[" << order << "] = " << rtn );
  return rtn;
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "SFE::internalFindC(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "findC() for eval point");
    evalPt()->findC(order, context);
  }

  if (order==1)
    {
      Deriv d(funcDeriv(this));
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }

  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_MSG3(verb,  tab << "SFE: C[" << order << "] = " << rtn );
  return rtn;
}


RCP<Array<Set<MultipleDeriv> > > SymbolicFuncElement
::internalDetermineR(const EvalContext& context,
                     const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab << "SFE::internalDetermineR() for "
                     << toString());
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "determineR() for eval point");
    evalPt()->determineR(context, RInput);

    SUNDANCE_MSG3(verb, tab1 << "SFE::internalDetermineR() for "
      << toString() << " delegating to EE");
    return EvaluatableExpr::internalDetermineR(context, RInput);
  }
}

bool SymbolicFuncElement::evalPtIsZero() const
{
  TEST_FOR_EXCEPTION(evalPt_.get()==0, InternalError,
                     "null evaluation point in SymbolicFuncElement::evalPtIsZero()");
  bool isZero = 0 != dynamic_cast<const ZeroExpr*>(evalPt());
  bool isTest = 0 != dynamic_cast<const TestFuncElement*>(this);
  return isZero || isTest;
}

void SymbolicFuncElement::substituteZero() const 
{
  evalPt_ = rcp(new ZeroExpr());
}

void SymbolicFuncElement
::substituteFunction(const RCP<DiscreteFuncElement>& u0) const
{
  evalPt_ = u0;
}



bool SymbolicFuncElement::isIndependentOf(const Expr& u) const 
{
  Expr uf = u.flatten();
  for (int i=0; i<uf.size(); i++)
  {
    const ExprBase* p = uf[i].ptr().get();
    const SymbolicFuncElement* f = dynamic_cast<const SymbolicFuncElement*>(p);
    TEST_FOR_EXCEPTION(f==0, InternalError, "expected a list of functions, "
      " got " << u);
    if (fid().dofID() == f->fid().dofID()) return false;
  }
  return true;
}


bool SymbolicFuncElement::isLinearForm(const Expr& u) const
{
  Expr uf = u.flatten();
  for (int i=0; i<uf.size(); i++)
  {
    const ExprBase* p = uf[i].ptr().get();
    const SymbolicFuncElement* f = dynamic_cast<const SymbolicFuncElement*>(p);
    TEST_FOR_EXCEPTION(f==0, InternalError, "expected a list of functions, "
      " got " << u);
    if (fid().dofID() == f->fid().dofID()) return true;
  }
  return false;
}


