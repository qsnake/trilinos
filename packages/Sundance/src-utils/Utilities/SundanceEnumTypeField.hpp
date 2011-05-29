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

#ifndef SUNDANCE_ENUMTYPEFIELD_HPP
#define SUNDANCE_ENUMTYPEFIELD_HPP

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"

namespace Sundance
{

template <typename T> class EnumTypeField
{
public:
  /** */
  EnumTypeField(const T& type) : type_(type) {}

  /** */
  void assertType(const T& reqType) const 
    {
      TEST_FOR_EXCEPTION(reqType != type(), RuntimeError, 
        "expected type=" << reqType << ", found type=" << type());
    }

  /** */
  void assertNotType(const T& tabooType) const 
    {
      TEST_FOR_EXCEPTION(tabooType == type(), RuntimeError, 
        "type=" << tabooType << " is unexpected in this context");
    }

  /** */
  const T& type() const {return type_;}

  /** */
  bool isType(const T& t) const {return type_==t;}
private:
  T type_;
};

}

#endif
