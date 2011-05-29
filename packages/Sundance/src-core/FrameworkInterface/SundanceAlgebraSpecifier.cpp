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

#include "SundanceAlgebraSpecifier.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace std;
using namespace Sundance;
using namespace Teuchos;




AlgebraSpecifier::AlgebraSpecifier()
  : EnumTypeField<AlgebraType>(ScalarAT), 
    direction_(-1)
{}

AlgebraSpecifier::AlgebraSpecifier(int direction)
  : EnumTypeField<AlgebraType>(CoordCompAT), 
    direction_(direction)
{}

AlgebraSpecifier::AlgebraSpecifier(
  const AlgebraType& at)
  : EnumTypeField<AlgebraType>(at), direction_(-1)
{
  assertNotType(CoordCompAT);
}

int AlgebraSpecifier::direction() const 
{
  assertType(CoordCompAT);
  return direction_;
}

bool AlgebraSpecifier::operator<(const AlgebraSpecifier& other) const 
{
  if (type() < other.type()) return true;
  if (type() > other.type()) return false;
  
  if (isCoordinateComponent())
    return direction() < other.direction();
  return false;
}

string AlgebraSpecifier::toString() const 
{
  TeuchosOStringStream os;
  os << *this;
  return os.str();
}

namespace std
{

ostream& operator<<(std::ostream& os, 
  const Sundance::AlgebraType& at)
{
  switch(at)
  {
    case ScalarAT:
      os << "Scalar";
      break;
    case VectorAT:
      os << "Vector";
      break;
    case CoordCompAT:
      os << "CoordComp";
      break;
    case NormalAT:
      os << "Normal";
      break;
    default:
      TEST_FOR_EXCEPT(1);
  }
  return os;
}

ostream& operator<<(std::ostream& os, const Sundance::AlgebraSpecifier& as)
{
  os << as.type();
  if (as.isCoordinateComponent()) os << "(d=" << as.direction() << ")";
  else os << "()";
  return os;
}

}

namespace Sundance
{

/** \relates AlgebraSpecifier */
AlgebraSpecifier vectorAlgebraSpec()
{
  return AlgebraSpecifier(VectorAT);
}

/** AlgebraSpecifier */
AlgebraSpecifier scalarAlgebraSpec()
{
  return AlgebraSpecifier(ScalarAT);
}

/** AlgebraSpecifier */
AlgebraSpecifier normalAlgebraSpec()
{
  return AlgebraSpecifier(NormalAT);
}

/** AlgebraSpecifier */
AlgebraSpecifier coordAlgebraSpec(int dir)
{
  return AlgebraSpecifier(dir);
}

}
