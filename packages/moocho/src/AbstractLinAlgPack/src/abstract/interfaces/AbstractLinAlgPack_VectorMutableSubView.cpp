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

#include "AbstractLinAlgPack_VectorMutableSubView.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

VectorMutableSubView::VectorMutableSubView( const vec_mut_ptr_t& vec, const Range1D& rng )
{
  this->initialize(vec,rng);
}

void VectorMutableSubView::initialize( const vec_mut_ptr_t& vec, const Range1D& rng )
{
  namespace rcp = MemMngPack;
  VectorSubView::initialize(vec,rng);
  full_vec_ = vec;
  this->has_changed();
}

void VectorMutableSubView::set_uninitialized()
{
  VectorSubView::set_uninitialized();
  full_vec_ = Teuchos::null;
  this->has_changed();
}

// Overriddend form Vector

Vector::vec_ptr_t VectorMutableSubView::sub_view( const Range1D& rng ) const
{
  return VectorSubView::sub_view(rng); // Had to override to resolve conflicit!
}

// Overridden from VectorMutable

void VectorMutableSubView::set_ele( index_type i, value_type val )
{
  space_impl().validate_range(Range1D(i,i));
  full_vec_->set_ele( space_impl().rng().lbound() + i - 1, val );
  this->has_changed();
}

VectorMutable::vec_mut_ptr_t
VectorMutableSubView::sub_view( const Range1D& rng_in )
{
  namespace rcp = MemMngPack;
  const size_type this_dim = this->dim();
  const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
  space_impl().validate_range(rng);
  if( rng.lbound() == 1 && rng.ubound() == this_dim )
    return Teuchos::rcp(this,false); // Do not own memory!
  // Below we will return an object that is not this entire object so
  // we must wipe out cache.
  this->has_changed();
  const index_type this_offset = space_impl().rng().lbound() - 1;
  return Teuchos::rcp(
    new VectorMutableSubView(
      full_vec_
      ,Range1D( 
        this_offset  + rng.lbound()
        ,this_offset + rng.ubound() )
      ) );
}

void VectorMutableSubView::get_sub_vector( const Range1D& rng_in, RTOpPack::MutableSubVector* sub_vec )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION( !sub_vec, std::logic_error ,"VectorMutableSubView::get_sub_vector(...): Error!" );
#endif
  const index_type this_dim = this->dim();
  const Range1D rng = RangePack::full_range(rng_in,1,this_dim);
  space_impl().validate_range(rng);
  const index_type this_offset = space_impl().rng().lbound() - 1;
  full_vec_->get_sub_vector( rng + this_offset, sub_vec );
  sub_vec->setGlobalOffset( sub_vec->globalOffset() - this_offset );
}

void VectorMutableSubView::commit_sub_vector( RTOpPack::MutableSubVector* sub_vec )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION( !sub_vec, std::logic_error, "VectorMutableSubView::commit_sub_vector(...): Error!" );
#endif
  const index_type this_offset = space_impl().rng().lbound() - 1;
  sub_vec->setGlobalOffset( sub_vec->globalOffset() + this_offset );
  full_vec_->commit_sub_vector( sub_vec );
  this->has_changed();
}

void VectorMutableSubView::set_sub_vector( const RTOpPack::SparseSubVector& sub_vec_in )
{
  const index_type            this_offset = space_impl().rng().lbound() - 1;
  RTOpPack::SparseSubVector   sub_vec = sub_vec_in;
  sub_vec.setGlobalOffset( sub_vec.globalOffset() + this_offset );
  full_vec_->set_sub_vector( sub_vec );
  this->has_changed();
}

} // end namespace AbstractLinAlgPack
