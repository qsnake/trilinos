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

#include <assert.h>

#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactorySerial.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "AbstractLinAlgPack_MultiVectorMutableDense.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef TEUCHOS_DEBUG
#define CLASS_MEMBER_PTRS \
const VectorSpaceSerial  *_this = this; \
const size_type *_dim = &dim_;
#else
#define CLASS_MEMBER_PTRS
#endif

namespace AbstractLinAlgPack {

VectorSpaceSerial::VectorSpaceSerial( size_type dim )
{
  CLASS_MEMBER_PTRS
  initialize(dim);
}

void VectorSpaceSerial::initialize( size_type dim )
{
  CLASS_MEMBER_PTRS
  dim_ = dim;
}

// Overridden from VectorSpace

bool VectorSpaceSerial::is_compatible(const VectorSpace& a_vec_space ) const
{
  CLASS_MEMBER_PTRS
  return this->dim() == a_vec_space.dim() && a_vec_space.is_in_core();
}

bool VectorSpaceSerial::is_in_core() const
{
  return true;
}

index_type VectorSpaceSerial::dim() const
{
  CLASS_MEMBER_PTRS
  return dim_;
}

VectorSpace::space_fcty_ptr_t
VectorSpaceSerial::small_vec_spc_fcty() const
{
  CLASS_MEMBER_PTRS
  return Teuchos::rcp(new VectorSpaceFactorySerial());
}

VectorSpace::space_ptr_t
VectorSpaceSerial::clone() const
{
  CLASS_MEMBER_PTRS
  namespace mmp = MemMngPack;
  return Teuchos::rcp( new VectorSpaceSerial( dim_	) );
}

VectorSpace::vec_mut_ptr_t
VectorSpaceSerial::create_member() const
{
  CLASS_MEMBER_PTRS
  namespace mmp = MemMngPack;
  return Teuchos::rcp(new VectorMutableDense(dim_));
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceSerial::create_members(size_type num_vecs) const
{
  CLASS_MEMBER_PTRS
  namespace mmp = MemMngPack;
  return Teuchos::rcp(new MultiVectorMutableDense(dim_,num_vecs));
}

VectorSpace::space_ptr_t
VectorSpaceSerial::sub_space(const Range1D& rng_in) const
{
  CLASS_MEMBER_PTRS
  namespace mmp = MemMngPack;
  const size_type this_dim = this->dim();
  const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    rng.ubound() > this_dim, std::out_of_range
    ,"VectorSpaceSerial::sub_view(...) : Error, "
    "rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1," << this_dim );
#endif
  if( rng == Range1D(1,this_dim) )
    return Teuchos::rcp( this, false );
  return Teuchos::rcp( new VectorSpaceSerial( rng.size() ) ); 
}

VectorSpace::space_ptr_t
VectorSpaceSerial::space(
  const GenPermMatrixSlice  &P
  ,BLAS_Cpp::Transp         P_trans
  ) const
{
  CLASS_MEMBER_PTRS
  return Teuchos::rcp( new VectorSpaceSerial( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ) ) ); 
}

} // end namespace AbstractLinAlgPack
