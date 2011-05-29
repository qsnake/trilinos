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

#ifndef SPARSE_LIN_ALG_PACK_SPARSE_BOUNDS_H
#define SPARSE_LIN_ALG_PACK_SPARSE_BOUNDS_H

#include "DenseLinAlgPack_DVectorClass.hpp"

namespace AbstractLinAlgPack {

/** \brief Iterate through a set of sparse bounds.
 *
 * Finish documentation.
 *
 * Allow default copy constructor and assignment operator.
 */
class sparse_bounds_itr {
private:
  enum EBound { LOWER, UPPER, BOTH };
public:

  /** \brief . */
  sparse_bounds_itr(	
    const DVectorSlice::const_iterator   &bl_begin
    ,const DVectorSlice::const_iterator  &bl_end
    ,const DVectorSlice::const_iterator  &bu_begin
    ,const DVectorSlice::const_iterator  &bu_end
    ,value_type                         big_bnd
    )
    :bl_itr_(bl_begin), bl_end_(bl_end), bu_itr_(bu_begin), bu_end_(bu_end)
    ,big_bnd_(big_bnd), index_(1)
  {
    if( !at_end() && ( *bl_itr_ <= -big_bnd_ && +big_bnd_ <= *bu_itr_ ) )
      this->operator++();
    else
      update();
  }
  /** \brief . */
  const value_type& big_bnd() const
  {    return big_bnd_; }
  /** \brief . */
  bool at_end() const
  {	return bl_itr_ == bl_end_; }
  /** \brief . */
  sparse_bounds_itr& operator++() {
    if(!at_end()) { ++bl_itr_; ++bu_itr_; ++index_; }
    for( ; !at_end() && ( *bl_itr_ <= -big_bnd_ && +big_bnd_ <= *bu_itr_ )
       ; ++bl_itr_, ++bu_itr_, ++index_ );
    update();
    return *this;
  }
  /** \brief . */
  index_type index() const
  {	return index_; }
  /** \brief . */
  value_type lbound() const
  {	return lbound_; }
  /** \brief . */
  value_type ubound() const
  {	return ubound_; }

private:
  DVectorSlice::const_iterator
    bl_itr_, bl_end_, bu_itr_, bu_end_;
  value_type
    big_bnd_, lbound_, ubound_;
  index_type
    index_;
  EBound
    at_bound_;

  void update() {
    if( bl_itr_ == bl_end_ ) {
      return;
    }
    else if( -big_bnd_ < *bl_itr_ && *bu_itr_ < +big_bnd_ ) {
      lbound_ = *bl_itr_;
      ubound_ = *bu_itr_;
      at_bound_ = BOTH;
    }
    else if( -big_bnd_ < *bl_itr_ ) {
      lbound_ = *bl_itr_;
      ubound_ = +big_bnd_;
      at_bound_ = LOWER;
    }
    else if( *bu_itr_ < +big_bnd_ ) {
      lbound_ = -big_bnd_;
      ubound_ = *bu_itr_;
      at_bound_ = UPPER;
    }
  }

  // not defined and not to be called
  sparse_bounds_itr();

};	// end class sparse_bounds_itr

}	// end namespace AbstractLinAlgPack

#endif // SPARSE_LIN_ALG_PACK_SPARSE_BOUNDS_H
