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

#ifndef SUNDANCE_MULTIINDEX_H
#define SUNDANCE_MULTIINDEX_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLObject.hpp"
#include <string>
#include <stdexcept>

namespace Sundance
{
using namespace Teuchos;

/**
 * An integer vector representing a multivariate derivative.
 */

class MultiIndex
{
public:
  /** constructs D(0,0,0) */
  MultiIndex();
  /** constructs a multiindex D(x,y,z) */
  MultiIndex(int x, int y, int z);

  /** */
  bool operator==(const MultiIndex& other) const ;

  /** */
  bool operator<(const MultiIndex& other) const ;

  /** */
  const int& operator[](int i) const {return m_[i];}

  /** */
  int& operator[](int i) {return m_[i];}

  /** */
  MultiIndex operator+(const MultiIndex& other) const ;

  /** */
  MultiIndex operator-(const MultiIndex& other) const ;

  /** */
  MultiIndex operator-() const ;

  /** */
  std::string toString() const ;

  /** */
  XMLObject toXML() const ;

  /** */
  int order() const ;

  /** */
  int firstOrderDirection() const ;

  /** */
  static int maxDim() {return 3;}

  /** */
  bool isValid() const ;

  /** */
  std::string coordForm() const ;
private:
  Array<int> m_;
};
}

namespace Teuchos
{

/** \relates Sundance::MultiIndex */
inline std::string toString(const Sundance::MultiIndex& h)
{return h.toString();}

}

namespace std
{
/** \relates Sundance::MultiIndex */
inline ostream& operator<<(std::ostream& os, 
  const Sundance::MultiIndex& h)
{
  os << h.toString();
  return os;
}
}

#endif
