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

#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorSpaceSubSpace.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_InnerProductDot.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

VectorSpace::VectorSpace( const inner_prod_ptr_t& inner_prod )
{
  this->inner_prod(inner_prod);
}

void VectorSpace::inner_prod( const inner_prod_ptr_t& inner_prod )
{
  if(inner_prod.get()) {
    inner_prod_ = inner_prod;
  } else {
    inner_prod_ = Teuchos::rcp(new InnerProductDot());
  }
}

const VectorSpace::inner_prod_ptr_t
VectorSpace::inner_prod() const
{
  return inner_prod_;
}

// Virtual functions with default implementations

bool VectorSpace::is_in_core() const
{
  return false;
}

VectorSpace::space_fcty_ptr_t
VectorSpace::small_vec_spc_fcty() const
{
  return Teuchos::null;
}

VectorSpace::vec_mut_ptr_t
VectorSpace::create_member(const value_type& alpha) const
{
  namespace mmp = MemMngPack;
  vec_mut_ptr_t vec = this->create_member();
  *vec = alpha;
  return vec;
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpace::create_members(size_type num_vecs) const
{
  return Teuchos::null;
}

VectorSpace::space_ptr_t
VectorSpace::sub_space(const Range1D& rng_in) const
{
  namespace mmp = MemMngPack;
  const index_type dim = this->dim();
  const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    rng.ubound() > dim, std::out_of_range
    ,"VectorSpace::sub_space(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
  if( rng.lbound() == 1 && rng.ubound() == dim )
    return space_ptr_t( this, false );
  return Teuchos::rcp(
    new VectorSpaceSubSpace(
      Teuchos::rcp( this, false )
      ,rng ) );
}

VectorSpace::space_ptr_t
VectorSpace::space(
  const GenPermMatrixSlice  &P
  ,BLAS_Cpp::Transp         P_trans
  ) const
{
  const index_type
    dim = BLAS_Cpp::rows( P.rows(), P.cols(), P_trans );
  space_fcty_ptr_t  vec_spc_fcty = this->small_vec_spc_fcty();
  if(vec_spc_fcty.get())
    return vec_spc_fcty->create_vec_spc(dim);
  return Teuchos::null;
}

// Overridden from AbstractFactory<>

VectorSpace::obj_ptr_t VectorSpace::create() const
{
  return this->create_member();
}

} // end namespace AbstractLinAlgPack
