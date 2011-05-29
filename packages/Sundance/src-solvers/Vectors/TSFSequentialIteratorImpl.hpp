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

#ifndef TSFSEQUENTIALITERATORIMPL_HPP
#define TSFSEQUENTIALITERATORIMPL_HPP


#include "TSFSequentialIteratorDecl.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "Thyra_SUNDIALS_Ops.hpp"
#include "TSFIndexableVector.hpp"


namespace TSFExtended
{

template <class Scalar> inline
bool SequentialIterator<Scalar>
::operator==(const SequentialIterator<Scalar>& other) const
{
  if (this->atEnd_ == other.atEnd_) return true;
  if (this->blockIndex_ != other.blockIndex_
    || this->globalIndex_ != other.globalIndex_
    || this->indexInCurrentBlock_ != other.indexInCurrentBlock_
    || *(this->space_) != *(other.space_)) return false;
  return true;
}



template <class Scalar> inline
SequentialIterator<Scalar> SequentialIterator<Scalar>::operator++(int)
{
  SequentialIterator<Scalar> old = *this;
  bool dataRemains = (this->space_)->advanceIndex(this->blockIndex_, 
    this->indexInCurrentBlock_, this->globalIndex_);
  if (!dataRemains) this->atEnd_ = true;
  return old;
}

template <class Scalar> inline
const VectorSpace<Scalar>& SequentialIterator<Scalar>::space() const
{
  return *(this->space_);
}


template <class Scalar> inline
SequentialIterator<Scalar>::SequentialIterator(
  const VectorSpace<Scalar>& space, 
  int blockIndex, int localIndex, int globalIndex
  )
  : space_(rcp(new VectorSpace<Scalar>(space))),
    blockIndex_(blockIndex),
    indexInCurrentBlock_(localIndex),
    globalIndex_(globalIndex),
    atEnd_(false)
{}

template <class Scalar> inline
SequentialIterator<Scalar>::SequentialIterator(
  const VectorSpace<Scalar>& space
  )
  : space_(rcp(new VectorSpace<Scalar>(space))),
    blockIndex_(-1),
    indexInCurrentBlock_(-1),
    globalIndex_(-1),
    atEnd_(true)
{}


template <class Scalar> inline
std::ostream& SequentialIterator<Scalar>::print(std::ostream& os) const 
{
  if (atEnd_) 
  {
    os << "[end]";
  }
  else
  {
    os << "[b=" << blockIndex_ << ", loc=" << indexInCurrentBlock_
       << ", glob=" << globalIndex_ 
       << "]";
  }
  return os;
}



}





#endif
