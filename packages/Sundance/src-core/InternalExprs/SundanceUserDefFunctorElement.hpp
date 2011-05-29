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

#ifndef SUNDANCE_USERDEFFUNCTORELEMENT_H
#define SUNDANCE_USERDEFFUNCTORELEMENT_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "Teuchos_Array.hpp"




namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/**
 * Scalar-valued element of a vector-valued functor
 */
class UserDefFunctorElement
{
public:
  /** ctor */
  UserDefFunctorElement(const RCP<const UserDefFunctor>& functor,
    int myIndex);

  /** */
  virtual ~UserDefFunctorElement(){;}

  /** */
  const std::string& name() const {return master_->name(myIndex());}

  /** */
  const std::string& masterName() const {return master_->name();}

  /** */
  void evalArgDerivs(int maxOrder, 
    const Array<double>& in,
    Array<double>& outDerivs) const ;

  /** */
  void getArgDerivIndices(const Array<int>& orders,
    Sundance::Map<MultiSet<int>, int>& varArgDerivs,
    Sundance::Map<MultiSet<int>, int>& constArgDerivs) const ;

  /** */
  int numArgs() const {return master_->domainDim();}

  /** */
  void reset() const {master_->reset();}

  /** Return the index of this element into the list-valued 
   * user defined op */
  int myIndex() const {return myIndex_;}

  /** */
  const UserDefFunctor* master() const 
    {return master_.get();}

  /** */
  int maxOrder() const {return master_->maxOrder();}


private:
  const RCP<const UserDefFunctor> master_;
  const int myIndex_;
};


}


#endif
