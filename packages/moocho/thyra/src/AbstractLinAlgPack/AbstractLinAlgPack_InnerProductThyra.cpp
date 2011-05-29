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

#include <stdexcept>

#include "AbstractLinAlgPack_InnerProductThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

InnerProductThyra::InnerProductThyra()
{}

InnerProductThyra::InnerProductThyra(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc
  )
{
  this->initialize(thyra_vec_spc);
}

void InnerProductThyra::initialize(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc
  )
{
  TEST_FOR_EXCEPTION(
    thyra_vec_spc.get()==NULL, std::invalid_argument
    ,"InnerProductThyra::initialize(thyra_vec_spc): Error!"
    );
  thyra_vec_spc_ = thyra_vec_spc;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > 
InnerProductThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > tmp_thyra_vec_spc = thyra_vec_spc_;
  thyra_vec_spc_ = Teuchos::null;
  return tmp_thyra_vec_spc;
}

// Overridden from InnerProduct

value_type InnerProductThyra::inner_prod(const Vector& v1, const Vector& v2) const
{
  using Teuchos::dyn_cast;
  return thyra_vec_spc_->scalarProd(
    *dyn_cast<const VectorMutableThyra>(v1).thyra_vec()
    ,*dyn_cast<const VectorMutableThyra>(v1).thyra_vec()
    );
}

} // end namespace AbstractLinAlgPack
