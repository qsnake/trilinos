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

#include "SundanceVectorBasisComponent.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

VectorBasisComponent::VectorBasisComponent(const BasisFamily& master, 
  int direction)
  : master_(master),
    direction_(direction)
{
  TEST_FOR_EXCEPTION(master_.ptr()->isScalarBasis(), RuntimeError,
    "scalar-valued basis " << master << " given as master argument to "
    "a vector basis component object. The master object must be a "
    "vector-valued basis");
}

bool VectorBasisComponent::lessThan(const BasisFamilyBase* other) const 
{
  TEST_FOR_EXCEPTION(
    (typeid(*this).before(typeid(*other)) 
      || typeid(*other).before(typeid(*this))),
    InternalError,
    "mismatched types: this=" << typeid(*this).name()
    << " and other=" << typeid(*other).name() 
    << " in BasisFamilyBase::lessThan(). This is most likely "
    "an internal bug, because the case of distinct types should have "
    "been dealt with before this point.");

  const VectorBasisComponent* p
    = dynamic_cast<const VectorBasisComponent*>(other);

  if (direction_ < p->direction_) return true;
  if (direction_ > p->direction_) return false;

  return (master_ < p->master_);
}
