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


#ifndef SUNDANCE_USERDEFOPEVALUATOR_H
#define SUNDANCE_USERDEFOPEVALUATOR_H

#include "SundanceDefs.hpp"

#include "SundanceUserDefFunctorElement.hpp"
#include "SundanceUserDefOpCommonEvaluator.hpp"
#include "SundanceChainRuleEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Sundance 
{

class UserDefOpElement;
class SymbolicFuncElementEvaluator;
class UserDefOpCommonEvaluator;
/**
 *
 */
class UserDefOpEvaluator : public ChainRuleEvaluator
{
public:
  /** */
  UserDefOpEvaluator(const UserDefOpElement* expr,
    const RCP<const UserDefOpCommonEvaluator>& commonEval,
    const EvalContext& context);

  /** */
  virtual ~UserDefOpEvaluator(){;}


  /** */
  virtual void 
  evalArgDerivs(const EvalManager& mgr,
    const Array<RCP<Array<double> > >& constArgRes,
    const Array<RCP<Array<RCP<EvalVector> > > >& vArgResults,
    Array<double>& constArgDerivs,
    Array<RCP<EvalVector> >& varArgDerivs) const ;

     
        
  /** */
  TEUCHOS_TIMER(evalTimer, "user defined nonlinear op evaluation");

  /** */
  void resetNumCalls() const ;

protected:

  Array<int> findRequiredOrders(const ExprWithChildren* expr, 
    const EvalContext& context) ;

  const UserDefFunctorElement* functor() const {return functor_;}

  const UserDefOpCommonEvaluator* commonEval() const 
    {return commonEval_.get();}

  int myIndex() const {return functor_->myIndex();}

private:
  Array<int> argValueIndex_;
  Array<int> argValueIsConstant_;
  const UserDefFunctorElement* functor_;
  RCP<const UserDefOpCommonEvaluator> commonEval_;
  int maxOrder_;
  int numVarArgDerivs_;
  int numConstArgDerivs_;
  bool allArgsAreConstant_;
}; 
}


#endif
