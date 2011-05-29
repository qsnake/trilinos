#if 0

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

#include "ConstrainedOptPack_MatrixHessianSuperBasicInitDiagonal.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Midynamic_cast_verbose.h"

namespace ConstrainedOptPack {

MatrixHessianSuperBasicInitDiagonal::MatrixHessianSuperBasicInitDiagonal()
  : B_RR_init_(NULL)
{}

void MatrixHessianSuperBasicInitDiagonal::initialize(
  size_type            n
  ,size_type           n_R
  ,const size_type     i_x_free[]
  ,const size_type     i_x_fixed[]
  ,const EBounds       bnd_fixed[]
  ,const B_RR_ptr_t&   B_RR_ptr
  ,const B_RX_ptr_t&   B_RX_ptr
  ,BLAS_Cpp::Transp    B_RX_trans
  ,const B_XX_ptr_t&   B_XX_ptr
  )
{
  using Teuchos::dyn_cast;

  // Validate the B_RR supports this interface
#ifdef _WINDOWS
  B_RR_init_ = &dynamic_cast<MatrixSymInitDiag&>(
    const_cast<MatrixSymWithOpFactorized&>(*B_RR_ptr)
    );
#else
  B_RR_init_ = &dyn_cast<MatrixSymInitDiag>(
    const_cast<MatrixSymWithOpFactorized&>(*B_RR_ptr)
    );
#endif

  MatrixHessianSuperBasic::initialize(
    n,n_R,i_x_free,i_x_fixed,bnd_fixed
    ,B_RR_ptr,B_RX_ptr,B_RX_trans,B_XX_ptr 
    );
}

// Overridden from MatrixSymInitDiag

void MatrixHessianSuperBasicInitDiagonal::init_identity(
  size_type n, value_type alpha )
{
  assert_initialized();
  B_RR_init_->init_identity(n,alpha);
  MatrixHessianSuperBasic::initialize(
    n,n,NULL,NULL,NULL
    ,this->B_RR_ptr()
    ,this->B_RX_ptr(),this->B_RX_trans()
    ,this->B_XX_ptr()
    );
}

void MatrixHessianSuperBasicInitDiagonal::init_diagonal(
  const DVectorSlice& diag )
{
  assert_initialized();
  B_RR_init_->init_diagonal(diag);
  MatrixHessianSuperBasic::initialize(
    diag.size(),diag.size(),NULL,NULL,NULL
    ,this->B_RR_ptr()
    ,this->B_RX_ptr(),this->B_RX_trans()
    ,this->B_XX_ptr()
    );
}

} // end namespace ConstrainedOptPack

#endif // 0
