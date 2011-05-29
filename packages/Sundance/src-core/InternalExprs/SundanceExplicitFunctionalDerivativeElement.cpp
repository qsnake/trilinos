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

#include "SundanceExplicitFunctionalDerivativeElement.hpp"

#include "SundanceEFDEEvaluator.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


ExplicitFunctionalDerivativeElement
::ExplicitFunctionalDerivativeElement(
  const RCP<ScalarExpr>& arg,
  const Deriv& fd
  )
  : UnaryExpr(arg), fd_(fd)
{
  Tabs tabs;
}

std::ostream& ExplicitFunctionalDerivativeElement
::toText(std::ostream& os, bool /* paren */) const 
{
	os << "FD[" << arg().toString() << ", " << fd_ << "]";
	return os;
}


XMLObject ExplicitFunctionalDerivativeElement
::toXML() const 
{
	XMLObject rtn("EFDE");
	rtn.addAttribute("df", fd_.toString());
  rtn.addChild(arg().toXML());

	return rtn;
}


Evaluator* ExplicitFunctionalDerivativeElement
::createEvaluator(const EvaluatableExpr* expr,
  const EvalContext& context) const
{
  return new EFDEEvaluator(dynamic_cast<const ExplicitFunctionalDerivativeElement*>(expr), context);
}




RCP<Array<Set<MultipleDeriv> > > 
ExplicitFunctionalDerivativeElement
::internalDetermineR(const EvalContext& context,
                           const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab0 << "ExplicitFunctionalDerivativeElement::internalDetermineR for=" << toString());
  SUNDANCE_MSG3(verb, tab0 << "RInput = " << RInput );

  RCP<Array<Set<MultipleDeriv> > > rtn 
    = rcp(new Array<Set<MultipleDeriv> >(RInput.size()));
  
  {
    Tabs tab1;
    for (int i=0; i<RInput.size(); i++)
      {
        Tabs tab2;
        const Set<MultipleDeriv>& Wi = findW(i, context);
        SUNDANCE_MSG5(verb,  tab2 << "W[" << i << "] = " << Wi );
        (*rtn)[i] = RInput[i].intersection(Wi);
      }

    Array<Set<MultipleDeriv> > RArg(RInput.size()+1);
    MultipleDeriv me(fd_);
    
    for (int order=1; order<=RInput.size(); order++)
      {
        Tabs tab2;
        SUNDANCE_MSG3(verb, tab2 << "order = " << order);
        if (RInput[order-1].size() == 0) continue;
        const Set<MultipleDeriv>& WArg = evaluatableArg()->findW(order, context);
        const Set<MultipleDeriv>& RMinus = (*rtn)[order-1];

        SUNDANCE_MSG3(verb, tab2 << "RInput = " << RInput[order-1]);
        SUNDANCE_MSG3(verb, tab2 << "FD times RInput = " 
          << setProduct(RMinus, makeSet(me)));
        
        RArg[order].merge(setProduct(RMinus, makeSet(me)).intersection(WArg));
      }
    SUNDANCE_MSG3(verb, tab1 << "RArg = " << RArg);
    
    SUNDANCE_MSG3(verb, tab1 << "calling determineR() for arg "
                       << evaluatableArg()->toString());
    evaluatableArg()->determineR(context, RArg);
  }
  SUNDANCE_MSG3(verb, tab0 << "R = " << (*rtn) );
  SUNDANCE_MSG3(verb, tab0 << "done with ExplicitFunctionalDerivativeElement::internalDetermineR for "
                     << toString());
  /* all done */  
  return rtn;
}


Set<MultipleDeriv> ExplicitFunctionalDerivativeElement
::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab0 
    << "ExplicitFunctionalDerivativeElement::internalFindW(order="
    << order << ") for " << toString());


  Set<MultipleDeriv> rtn ;

  if (order < 3)
  {
    Tabs tab1;
    const Set<MultipleDeriv>& WArgPlus 
      = evaluatableArg()->findW(order+1, context);
    
    SUNDANCE_MSG5(verb, tab1 << "WArgPlus = " << WArgPlus);
    MultipleDeriv me(fd_);
    Set<MultipleDeriv> WargPlusOslashFD = setDivision(WArgPlus, makeSet(me));
    SUNDANCE_MSG5(verb, tab1 << "WArgPlus / fd = " 
      << WargPlusOslashFD);
    rtn = WargPlusOslashFD;
  }
  SUNDANCE_MSG3(verb, tab0 << "W[" << order << "]=" << rtn);
  SUNDANCE_MSG3(verb, tab0 << "done with ExplicitFunctionalDerivativeElement::internalFindW(" << order << ") for "
                     << toString());

  return rtn;
}


Set<MultipleDeriv> ExplicitFunctionalDerivativeElement
::internalFindV(int order, const EvalContext& context) const
{
  Tabs tabs;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tabs << "ExplicitFunctionalDerivativeElement::internalFindV(" << order << ") for " 
                     << toString());

  Set<MultipleDeriv> rtn;
  
  {
    Tabs tab1;
    SUNDANCE_MSG5(verb, tab1 << "finding R");
    const Set<MultipleDeriv>& R = findR(order, context);
    SUNDANCE_MSG5(verb, tab1 << "finding C");
    const Set<MultipleDeriv>& C = findC(order, context);
    rtn = R.setDifference(C);
  }
  SUNDANCE_MSG3(verb, tabs << "V[" << order << "]=" << rtn);
  SUNDANCE_MSG3(verb, tabs << "done with ExplicitFunctionalDerivativeElement::internalFindV(" << order << ") for "
    << toString());
  
  return rtn;
}


Set<MultipleDeriv> ExplicitFunctionalDerivativeElement
::internalFindC(int order, const EvalContext& context) const
{
  Tabs tabs;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tabs << "ExplicitFunctionalDerivativeElement::internalFindC() for " 
    << toString());
  Set<MultipleDeriv> rtn ;
  if (order < 3) 
  {
    Tabs tab1;
    SUNDANCE_MSG5(verb, tab1 << "finding R");
    const Set<MultipleDeriv>& R = findR(order, context);
    SUNDANCE_MSG5(verb, tab1 << "R=" << R);

    SUNDANCE_MSG5(verb, tab1 << "finding C for arg");
    const Set<MultipleDeriv>& argC 
      = evaluatableArg()->findC(order+1, context);
    SUNDANCE_MSG5(verb, tab1 << "argC=" << argC);

    MultipleDeriv me(fd_);
    Set<MultipleDeriv> tmp = setDivision(argC, makeSet(me));
    rtn = tmp.intersection(R);
  }

  SUNDANCE_MSG3(verb, tabs << "C[" << order << "]=" << rtn);
  SUNDANCE_MSG3(verb, tabs << "done with ExplicitFunctionalDerivativeElement::internalFindC for "
    << toString());
  return rtn;
}




bool ExplicitFunctionalDerivativeElement
::lessThan(const ScalarExpr* other) const
{
  const ExplicitFunctionalDerivativeElement* e 
    = dynamic_cast<const ExplicitFunctionalDerivativeElement*>(other);
  TEST_FOR_EXCEPTION(e==0, InternalError, "cast should never fail at this point");
  
  if (fd_ < e->fd_) return true;
  if (e->fd_ < fd_) return false;
  
  return ExprWithChildren::lessThan(other);
}

