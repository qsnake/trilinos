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

#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace AbstractLinAlgPack {

BasisSystem::BasisSystem(
  const mat_sym_fcty_ptr_t             &factory_transDtD
  ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
  )
{
  this->initialize(factory_transDtD,factory_S);
}

void BasisSystem::initialize(
  const mat_sym_fcty_ptr_t             &factory_transDtD
  ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
  )
{
  factory_transDtD_ = factory_transDtD;
  factory_S_        = factory_S;
}

Range1D BasisSystem::equ_decomp() const
{
  const size_type r = this->var_dep().size();
  return r ? Range1D(1,r) : Range1D::Invalid;
}

Range1D BasisSystem::equ_undecomp() const
{
  return Range1D::Invalid;
}

const BasisSystem::mat_fcty_ptr_t BasisSystem::factory_GcUP() const
{
  return Teuchos::null;
}

const BasisSystem::mat_sym_fcty_ptr_t
BasisSystem::factory_transDtD() const
{
  return factory_transDtD_;
}
  
const BasisSystem::mat_sym_nonsing_fcty_ptr_t
BasisSystem::factory_S() const
{
  return factory_S_;
}

} // end namespace AbstractLinAlgPack
