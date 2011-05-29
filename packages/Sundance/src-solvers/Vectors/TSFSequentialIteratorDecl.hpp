/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
 /* @HEADER@ */

#ifndef TSFSEQUENTIALITERATORDECL_HPP
#define TSFSEQUENTIALITERATORDECL_HPP

#include "SundanceDefs.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

namespace TSFExtended
{

using Teuchos::RCP;

/* Forward declare VectorSpace */
template <class Scalar> 
class VectorSpace;

/* Forward declare Vector */
template <class Scalar> 
class Vector;


/** */
template <class Scalar>
class SequentialIterator
{
public:

  friend class VectorSpace<Scalar>;
  friend class Vector<Scalar>;

  /** */
  SequentialIterator() : space_(), blockIndex_(0), indexInCurrentBlock_(0),
                         atEnd_(true){}

  /** */
  bool operator==(const SequentialIterator<Scalar>& other) const ;

  /** */
  bool operator!=(const SequentialIterator<Scalar>& other) const 
    {
      return !operator==(other);
    }

  /** */
  SequentialIterator<Scalar> operator++(int);

  /** */
  const VectorSpace<Scalar>& space() const ;

  
  /** */
  std::ostream& print(std::ostream& os) const ;

  /** */
  const OrdType& globalIndex() const {return globalIndex_;}

private:

  /** */
  const OrdType& indexInBlock() const {return indexInCurrentBlock_;}

  /** */
  const OrdType& blockIndex() const {return blockIndex_;}

  /** Constructor is private: the construction is always done inside
   * the begin and end methods of vector space. */
  SequentialIterator(const VectorSpace<Scalar>& space, 
    int blockIndex, int indexInCurrentBlock, int globalIndex);

  /** Constructor is private: the construction is always done inside
   * the begin and end methods of vector space. */
  SequentialIterator(const VectorSpace<Scalar>& space);

  /* This is odd but necessary: we store the VectorSpace handle in an RCP
   * because we can only forward declare VectorSpace at this point. */
  RCP<VectorSpace<Scalar> > space_;
  OrdType blockIndex_;
  OrdType indexInCurrentBlock_;
  OrdType globalIndex_;
  bool atEnd_;
};

}

namespace TSFExtended
{
template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, 
  const TSFExtended::SequentialIterator<Scalar>& i)
{
  return i.print(os);
}
}


#endif
