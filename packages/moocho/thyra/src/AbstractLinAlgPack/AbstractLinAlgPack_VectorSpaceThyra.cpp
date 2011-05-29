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

#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactoryThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MultiVectorMutableThyra.hpp"
#include "AbstractLinAlgPack_InnerProductThyra.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

VectorSpaceThyra::VectorSpaceThyra()
{}

VectorSpaceThyra::VectorSpaceThyra(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
  ,const inner_prod_ptr_t                                                  &inner_prod
  )
{
  this->initialize(thyra_vec_spc,inner_prod);
}

void VectorSpaceThyra::initialize(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc,
  const inner_prod_ptr_t &inner_prod
  )
{
  namespace mmp = MemMngPack;
  TEST_FOR_EXCEPTION(
    thyra_vec_spc.get()==NULL, std::invalid_argument
    ,"VectorSpaceThyra::initialize(thyra_vec_spc): Error!"
    );
  thyra_vec_spc_ = thyra_vec_spc;
  if(inner_prod.get())
    this->inner_prod(inner_prod);
  else
    this->inner_prod(Teuchos::rcp(new InnerProductThyra(thyra_vec_spc)));
}

Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >
VectorSpaceThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > tmp_thyra_vec_spc = thyra_vec_spc_;
  thyra_vec_spc_ = Teuchos::null;
  return tmp_thyra_vec_spc;
}

// Overridden from VectorSpace

VectorSpace::space_ptr_t
VectorSpaceThyra::clone() const
{
  return Teuchos::rcp(new VectorSpaceThyra(thyra_vec_spc_->clone()));
}

bool VectorSpaceThyra::is_compatible(const VectorSpace& vec_spc ) const
{
  if( this->dim()==vec_spc.dim() && this->is_in_core() && vec_spc.is_in_core() )
    return true;
  const VectorSpaceThyra
    *thyra_vec_spc = dynamic_cast<const VectorSpaceThyra*>(&vec_spc);
  if( thyra_vec_spc->thyra_vec_spc()->isCompatible(*thyra_vec_spc_) )
    return true;
  return false;
}

bool VectorSpaceThyra::is_in_core() const
{
  return thyra_vec_spc_->hasInCoreView();
}

index_type VectorSpaceThyra::dim() const
{
  return thyra_vec_spc_->dim();
}

VectorSpace::vec_mut_ptr_t
VectorSpaceThyra::create_member() const
{
  return Teuchos::rcp(new VectorMutableThyra(Thyra::createMember(thyra_vec_spc_)));
}

VectorSpace::space_fcty_ptr_t
VectorSpaceThyra::small_vec_spc_fcty() const
{
  return Teuchos::rcp(new VectorSpaceFactoryThyra(thyra_vec_spc_->smallVecSpcFcty()));
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceThyra::create_members(size_type num_vecs) const
{
  return Teuchos::rcp(new MultiVectorMutableThyra(Thyra::createMembers(thyra_vec_spc_,num_vecs)));
}

} // end namespace AbstractLinAlgPack
