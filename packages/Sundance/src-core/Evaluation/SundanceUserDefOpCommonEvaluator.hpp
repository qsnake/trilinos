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


#ifndef SUNDANCE_USERDEFOPCOMMONEVALUATOR_H
#define SUNDANCE_USERDEFOPCOMMONEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalVector.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Sundance {

class UserDefOpElement;
class EvalContext;

/** 
 * UserDefOpCommonEvaluator provides a single evaluation point for all
 * components of a vector-valued functor. 
 */
class UserDefOpCommonEvaluator : public ObjectWithClassVerbosity<Evaluator>
{
public:
  /** */
  UserDefOpCommonEvaluator(const UserDefFunctor* op,
    const UserDefOpElement* expr,
    const EvalContext& context);

  /** */
  virtual ~UserDefOpCommonEvaluator(){;}

  /** Evaluate all vector components at the specified argument values. */
  void evalAllComponents(const EvalManager& mgr,
    const Array<RCP<Array<double> > >& constArgDerivVals,
    const Array<RCP<Array<RCP<EvalVector> > > >& vArgDerivVals) const ;

      
      
  /** Get the cached constant-valued argument derivative values
   * for the specified vector component */
  const Array<double>& constArgDerivCache(int elemIndex) const 
    {return *(constArgDerivCache_[elemIndex]);}

  /** Get the cached vector-valued argument derivative values
   * for the specified vector component */
  const Array<RCP<EvalVector> >& varArgDerivCache(int elemIndex) const 
    {return varArgDerivCache_[elemIndex];}

  /** Indicate whether the cached argument derivative values are valid */
  bool cacheIsValid() const {return cacheIsValid_;}

  /** Mark the cache as valid.  */
  void markCacheAsValid() const {cacheIsValid_ = true;}

  /** Mark the cache as invalid */
  void markCacheAsInvalid() const {cacheIsValid_ = false;}

  /** */
  int maxDiffOrder() const {return maxOrder_;}

  /** */
  void updateMaxOrder(int maxOrder) const {maxOrder_ = maxOrder;}


protected:

  /** Access the functor that implements this operator */
  const UserDefFunctor* functor() const {return functor_;}

private:
  mutable int maxOrder_;
  Array<int> argValueIndex_;
  Array<int> argValueIsConstant_;
  mutable Array<RCP<Array<double> > > constArgDerivCache_;
  mutable Array<Array<RCP<EvalVector> > > varArgDerivCache_;
  mutable bool cacheIsValid_;
  const UserDefFunctor* functor_;
      
}; 
}


#endif
