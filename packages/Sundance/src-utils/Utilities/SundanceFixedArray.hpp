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

#ifndef SUNDANCE_FIXEDARRAY_H
#define SUNDANCE_FIXEDARRAY_H


#include "SundanceDefs.hpp"
#include <iostream>
#include <algorithm>

namespace Sundance
{

/**
 * A simple fixed-size array 
 */
template <int N, class T> class FixedArray
{
public:
  /** */
  FixedArray(){;}

  /** */
  const T& operator[](int i) const 
    {
      TEST_FOR_EXCEPT(i<0);
      TEST_FOR_EXCEPT(i>=N);
      return data_[i];
    }

  /** */
  T& operator[](int i) 
    {
      TEST_FOR_EXCEPT(i<0);
      TEST_FOR_EXCEPT(i>=N);
      return data_[i];
    }

  /** */
  int size() const {return N;}

  /** */
  const T* begin() const {return &(data_[0]);}

  /** */
  const T* end() const {return begin()+N;}

  /** */
  T* begin() {return data_;}
  
  /** */
  T* end() {return data_+N;}

  /** */
  bool operator<(const FixedArray<N, T>& other) const
    {
      return std::lexicographical_compare(begin(), end(), 
        other.begin(), other.end());
    }

private:
  T data_[N];
};



/** \relate FixedArray */
template<class T>
FixedArray<2,T> makeFA(const T& a, const T& b)
{
  FixedArray<2,T> rtn;
  rtn[0] = a;
  rtn[1] = b;
  return rtn;
}

/** \relate FixedArray */
template<class T>
FixedArray<3,T> makeFA(const T& a, const T& b, const T& c)
{
  FixedArray<3,T> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  return rtn;
}


/** \relate FixedArray */
template<class T>
FixedArray<4,T> makeFA(const T& a, const T& b, const T& c, const T& d)
{
  FixedArray<4,T> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[2] = d;
  return rtn;
}

}


namespace std
{
/** \relates FixedArray */
template <int N, class T> 
ostream& operator<<(std::ostream& os, const Sundance::FixedArray<N, T>& f)
{
  os << "{";
  for (int i=0; i<N; i++) 
  {
    if (i>0) os << ", ";
    os << f[i];
  }
  os << "}";
  return os;
}

}

#endif
