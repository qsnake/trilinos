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

#include <ostream>

#include "AbstractLinAlgPack_MatrixSymIdent.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace AbstractLinAlgPack {

// Constructors/initalizers

MatrixSymIdent::MatrixSymIdent(
  const VectorSpace::space_ptr_t&          vec_space
  ,const value_type                        scale
  )
{
  this->initialize(vec_space,scale);
}

void MatrixSymIdent::initialize(
  const VectorSpace::space_ptr_t&          vec_space
  ,const value_type                        scale
  )
{
  vec_space_ = vec_space;
  scale_     = scale;
}

// Overridden from MatrixBase

size_type MatrixSymIdent::rows() const
{
  return vec_space_.get() ? vec_space_->dim() : 0;
}

size_type MatrixSymIdent::nz() const
{
  return vec_space_.get() ? vec_space_->dim() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixSymIdent::space_cols() const {
  return *vec_space_;
}

std::ostream& MatrixSymIdent::output(std::ostream& out) const
{
  out << "Identity matrix of dimension " << rows() << " x " << rows() << std::endl;
  return out;
}

void MatrixSymIdent::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const Vector& x, value_type b
  ) const
{
  AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, BLAS_Cpp::no_trans, x );
  Vt_S(y,b);
    Vp_StV(y,a*scale_,x);
}

// Overridden from MatrixNonsing

void MatrixSymIdent::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp M_trans, const Vector& x
  ) const
{
  AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, BLAS_Cpp::no_trans, x );
  LinAlgOpPack::V_StV(y,1.0/scale_,x);
}

} // end namespace AbstractLinAlgPack
