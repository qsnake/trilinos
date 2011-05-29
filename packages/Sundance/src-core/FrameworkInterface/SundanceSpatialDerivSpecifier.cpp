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

#include "SundanceSpatialDerivSpecifier.hpp"

using namespace Sundance;
using namespace Sundance;



SpatialDerivSpecifier::SpatialDerivSpecifier()
  : EnumTypeField<SpatialDerivType>(IdentitySDT), 
    mi_(), normalDerivOrder_(-1)
{}

SpatialDerivSpecifier::SpatialDerivSpecifier(const MultiIndex& mi)
  : EnumTypeField<SpatialDerivType>(PartialSDT), 
    mi_(mi), normalDerivOrder_(-1)
{}

SpatialDerivSpecifier::SpatialDerivSpecifier(
  const SpatialDerivType& sdt,
  int order)
  : EnumTypeField<SpatialDerivType>(sdt),
    mi_(), normalDerivOrder_(order)
{
  assertNotType(PartialSDT);

  if (order > 0)
  {
    assertNotType(DivSDT);
  }
}

const MultiIndex& SpatialDerivSpecifier::mi() const
{
  assertNotType(DivSDT);
  assertNotType(NormalSDT);
  return mi_;
}

bool SpatialDerivSpecifier::isDivergence() const
{
  return isType(DivSDT);
}

bool SpatialDerivSpecifier::isNormal() const
{
  return isType(NormalSDT);
}

bool SpatialDerivSpecifier::isPartial() const
{
  return isType(PartialSDT);
}

bool SpatialDerivSpecifier::isIdentity() const
{
  return isType(IdentitySDT)
    || (isPartial() && mi().order()==0)
    || (isNormal() && normalDerivOrder()==0);
}

int SpatialDerivSpecifier::normalDerivOrder() const
{
  assertType(NormalSDT);
  return normalDerivOrder_;
}

int SpatialDerivSpecifier::derivOrder() const
{
  if (isDivergence()) return 1;
  if (isPartial()) return mi_.order();
  if (isNormal()) return normalDerivOrder_;
  return 0;
}


std::string SpatialDerivSpecifier::toString() const 
{
  TeuchosOStringStream os;
  os << *this;
  return os.str();
}



bool SpatialDerivSpecifier::operator<(const SpatialDerivSpecifier& other) const
{
  if (type() < other.type()) return true;
  if (type() > other.type()) return false;

  if (isPartial()) return mi() < other.mi();
  if (isNormal()) return normalDerivOrder() < other.normalDerivOrder();

  return false;
}


SpatialDerivSpecifier SpatialDerivSpecifier::derivWrtMultiIndex(const MultiIndex& mi) const 
{
  if (isPartial() || isIdentity()) 
  {
    return SpatialDerivSpecifier(mi_+mi);
  }
  else if (mi.order()==0)
  {
    return *this;
  }
  else
  {
    TEST_FOR_EXCEPTION(true, InternalError, "cannot take an arbitrary "
      "spatial derivative of SDS=" << *this);
    return *this; // -Wall
  }
  
}


namespace std
{
/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivSpecifier& sds)
{
  os << sds.type();
  if (sds.isPartial()) os << "(d=" << sds.mi() << ")";
  else os << "()";
  return os;
}


/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivType& sdt)
{
  static Array<string> names = tuple<string>("Identity", "Partial", "Normal", "Divergence");
  os << names[sdt];
  return os;
}

}
