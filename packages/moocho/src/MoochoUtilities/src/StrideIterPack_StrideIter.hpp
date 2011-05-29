// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STRIDE_ITER_H
#define STRIDE_ITER_H

#include "Moocho_ConfigDefs.hpp"

namespace StrideIterPack {

/** \brief C++ Standard Library compatable iterator class for accesing nonunit stride arrays of data.
 *
 * Random access iterator adaptor for non-unit strides.
 * This iterator does not have any range checking for maximum efficency.  It was desinged
 * with the purpose of allowing a slice from an array (BLAS matrix row) to be used.  It,
 * however can be used with any random access iterator.  There are several
 * \ref stride_iter_funcs_grp "non-member functions" that can be used with this class.
 */
template<class T_iterator_type, class T_value_type, class T_reference_type
  , class T_pointer_type, class T_difference_type>
class stride_iter {
public:
  /** @name typedefs */
  //@{

  /** \brief . */
  typedef	std::random_access_iterator_tag						iterator_category;
  /** \brief . */
  typedef	T_iterator_type										iterator_type;
  /** \brief . */
  typedef	T_value_type										value_type;
  /** \brief . */
  typedef T_reference_type									reference;
  /** \brief . */
  typedef T_pointer_type										pointer;
  /** \brief . */
  typedef	T_difference_type									difference_type;
  
  //@}
    
  /** @name Constructors.  Uses default copy constructor and assignment operator */
  //@{

  /// constructs to a null iterator
  stride_iter() : current_(0), stride_(0)
  {}
  /// constructs to the standard iterator (behaves as) (increment 1).  Allows conversion from iterator.
  stride_iter(iterator_type current) :  current_(current), stride_(1)
  {}	
  /// constructs to the desired sliced iterator
  stride_iter(iterator_type current, difference_type stride) :  current_(current), stride_(stride)
  {}
  /// convert type of iterators (mainly for non-const to const)
  template<class Iter, class Val, class Ref, class Ptr, class Diff>
  stride_iter(const stride_iter<Iter,Val,Ref,Ptr,Diff>& rhs) :current_(rhs.current()),stride_(rhs.stride())
  {}
  /// assign different types of iterators (mainly for non-const to const)
  template<class Iter, class Val, class Ref, class Ptr, class Diff>
  stride_iter & operator=(const stride_iter<Iter,Val,Ref,Ptr,Diff>& rhs) {
    current_ = rhs.current(); stride_ = rhs.stride(); return *this;
  }

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  reference		operator*()	const {
    return *current_;
  }
  /** \brief . */
  pointer		operator->() const {
    return  current_;
  }
  /** \brief . */
  reference		operator[](difference_type n)	const {
    return current_[n * stride_];
  }

  //@}

  /** @name Incrementation */
  //@{

  /// ++itr
  stride_iter&		operator++() {
    current_ += stride_;
    return *this;
  }
  /// itr++
  const stride_iter	operator++(int) {
    stride_iter tmp = *this;
    ++*this; return tmp;
  }
  /// --itr
  stride_iter&		operator--() {
    current_ -= stride_;
    return *this;
  }
  /// itr--
  const stride_iter	operator--(int) {
    stride_iter tmp = *this;
    --*this; return tmp;
  }
  /// itr + n 
  stride_iter			operator+(difference_type n) {
    return stride_iter(current_ + n * stride_, stride_);
  }
  /// itr + n
  const stride_iter	operator+(difference_type n) const {
    return stride_iter(current_ + n * stride_, stride_);
  }
  /// itr += n
  stride_iter&		operator+=(difference_type n) {
    current_ += n * stride_;
    return *this;
  }
  /// itr - n
  stride_iter			operator-(difference_type n) {
    return stride_iter(current_ - n * stride_, stride_);
  }
  /// itr - n
  const stride_iter	operator-(difference_type n) const {
    return stride_iter(current_ - n * stride_, stride_);
  }
  /// itr -= n
  stride_iter&		operator-=(difference_type n) {
    current_ -= n * stride_;
    return *this;
  }

  //@}

  /// distance = itr1 - itr2 (distance between elements)
  difference_type		operator-(const stride_iter& itr) const	{
    return (current_ - itr.current_)/stride_;
  }

  /** @name Implimentation access */
  //@{

  /** \brief . */
  iterator_type	current() const {
    return current_;
  }

  /** \brief . */
  difference_type	stride() const {
    return stride_;
  }

  //@}

private:
  /** \brief . */
  iterator_type			current_;
  /** \brief . */
  difference_type			stride_;
};

// RAB: 2001/11/16: For some reason this \defgroup causing doxygen to crash
// when compiling RTOpPack?

/* \defgroup stride_iter_funcs_grp Non-member functions for StrideIterPack::stride_itr
 * \ingroup Misc_grp
 */
// @{

/// Allow difference_type as lhs argument in n + itr
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline const stride_iter<Iter,Val,Ref,Ptr,Diff> operator+(Diff n
  , const stride_iter<Iter,Val,Ref,Ptr,Diff> itr)	
{	
  return itr + n;
}

/// itr1 < itr2
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline bool operator<(const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr1, 
            const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr2)
{	
  return itr1.operator->() < itr2.operator->();
}
/// itr1 <= itr2
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline bool operator<=(const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr1, 
             const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr2)
{
  return itr1.operator->() <= itr2.operator->();
}
/// itr1 > itr 2
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline bool operator>(const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr1, 
            const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr2)
{
  return itr1.operator->() > itr2.operator->();
}
/// itr1 >= itr2
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline bool operator>=(const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr1, 
             const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr2)
{
  return itr1.operator->() >= itr2.operator->();
}
/// itr1 == itr2
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline bool operator==(const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr1, 
             const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr2)
{
  return itr1.operator->() == itr2.operator->();
}
/// itr1 != itr2
template<class Iter, class Val, class Ref, class Ptr, class Diff>
inline bool operator!=(const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr1, 
             const stride_iter<Iter,Val,Ref,Ptr,Diff>& itr2)
{
  return itr1.operator->() != itr2.operator->();
}

// @}

}	// end namespace StrideIterPack

#endif // STRIDE_ITER_H
