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

#include "SundanceFunctionIdentifier.hpp"
#include "SundanceOrderedTuple.hpp"

using namespace Sundance;
using namespace Sundance;

FunctionIdentifier::FunctionIdentifier()
  : dofID_(-1), algSpec_(ScalarAT)
{}

FunctionIdentifier::FunctionIdentifier(
  const AlgebraSpecifier& algSpec)
  : dofID_(nextID()), algSpec_(algSpec)
{}

FunctionIdentifier::FunctionIdentifier(
  const FunctionIdentifier* parent,
  const AlgebraSpecifier& algSpec)
  : dofID_(-1), algSpec_(algSpec)
{
  /* make sure the parent exists */
  TEST_FOR_EXCEPT(parent==0);
  dofID_ = parent->dofID();

  /* check for various stupid cases that should never happen */
  TEST_FOR_EXCEPTION(!parent->algSpec_.isVector(), InternalError,
    "attempted to form a function ID for a component of a non-vector object:"
    "parent=" << parent->toString() << " component spec=" << algSpec);
  TEST_FOR_EXCEPTION(algSpec.isVector() || algSpec.isScalar(), RuntimeError,
    "attempted to define a vector or scalar as a component of another object."
    "parent=" <<  parent->toString() << " component spec=" << algSpec);
}



int FunctionIdentifier::componentIndex() const 
{
  TEST_FOR_EXCEPTION(!(algSpec_.isCoordinateComponent() || algSpec_.isScalar()), InternalError,
    "attempted to find component index for a FID that is not a "
    "scalar or a coordinate component of a vector");
  if (algSpec_.isScalar()) return 0;
  return algSpec_.direction();
}

string FunctionIdentifier::toString() const 
{
  TeuchosOStringStream os;
  os << *this;
  return os.str();
}

FunctionIdentifier FunctionIdentifier::createNormal() const
{
  TEST_FOR_EXCEPTION(!isVector(), InternalError,
    "attempted to find normal component of a FID that is not a vector");

  return FunctionIdentifier(this, normalAlgebraSpec());
}

bool FunctionIdentifier::operator<(const FunctionIdentifier& other) const 
{
  OrderedPair<int, AlgebraSpecifier> me(dofID_, algSpec_);
  OrderedPair<int, AlgebraSpecifier> you(other.dofID_, other.algSpec_);
  return me < you;
}

FunctionIdentifier FunctionIdentifier::createComponent(int d) const
{
  TEST_FOR_EXCEPTION(!isVector(), InternalError,
    "attempted to find component of a FID that is not a vector");

  return FunctionIdentifier(this, coordAlgebraSpec(d));
}

namespace Sundance
{

FunctionIdentifier makeFuncID(int tensorOrder)
{
  if (tensorOrder==0)
    return FunctionIdentifier(scalarAlgebraSpec());
  else if (tensorOrder==1)
    return FunctionIdentifier(vectorAlgebraSpec());
  else
    TEST_FOR_EXCEPT(true);
  return FunctionIdentifier(scalarAlgebraSpec()); // -Wall
}

}


namespace std
{

ostream& operator<<(std::ostream& os, 
  const Sundance::FunctionIdentifier& fid)
{
  os << "FuncID(dofID=" << fid.dofID() << ", component type=" << fid.algSpec() << ")";
  return os;
}

}

