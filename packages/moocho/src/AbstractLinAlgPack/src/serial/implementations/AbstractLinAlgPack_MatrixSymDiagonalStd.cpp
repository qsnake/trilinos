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

#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace AbstractLinAlgPack {

MatrixSymDiagStd::MatrixSymDiagStd()
{}

DVectorSlice MatrixSymDiagStd::diag()
{
  return diag_();
}

const DVectorSlice MatrixSymDiagStd::diag() const
{
  return diag_();
}

// Overridden from MatrixSymInitDiag

void MatrixSymDiagStd::init_identity( size_type n, value_type alpha = 1.0 )
{
  diag_.resize(n);
  if(n)
    diag_ = alpha;
}

void MatrixSymDiagStd::init_diagonal( const DVectorSlice& diag )
{
  diag_ = diag;
}

// Overridden from Matrix

size_type MatrixSymDiagStd::rows() const
{
  return diag_.size();
}

// Overridden from MatrixOp

void MatrixSymDiagStd::Mp_StM(
  DMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const
{
  using DenseLinAlgPack::Vp_StV;
  // Assert that the dimensions match up.
  DenseLinAlgPack::Mp_M_assert_sizes( gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
    , rows(), cols(), trans_rhs );
  // Add to the diagonal
  Vp_StV( &gms_lhs->diag(), alpha, diag_ );
}

void MatrixSymDiagStd::Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2, value_type beta) const
{
  const DVectorSlice diag = this->diag();
  size_type n = diag.size();

  //
  // y = b*y + a * op(A) * x
  //
  DenseLinAlgPack::Vp_MtV_assert_sizes(
    vs_lhs->size(), n, n, trans_rhs1, vs_rhs2.size() );
  //
  // A is symmetric and diagonal A = diag(diag) so:
  //
  // y(j) += a * diag(j) * x(j), for j = 1...n
  //
  if( vs_rhs2.stride() == 1 && vs_lhs->stride() == 1 ) {
    // Optimized implementation
    const value_type
      *d_itr      = diag.raw_ptr(),
      *x_itr      = vs_rhs2.raw_ptr();
    value_type
      *y_itr      = vs_lhs->raw_ptr(),
      *y_end      = y_itr + vs_lhs->size();

    if( beta == 0.0 ) {
      while( y_itr != y_end )
        *y_itr++ = alpha * (*d_itr++) * (*x_itr++);
    }
    else if( beta == 1.0 ) {
      while( y_itr != y_end )
        *y_itr++ += alpha * (*d_itr++) * (*x_itr++);
    }
    else {
      for( ; y_itr != y_end; ++y_itr )
        *y_itr = beta * (*y_itr) + alpha * (*d_itr++) * (*x_itr++);
    }
  }
  else {
    // Generic implementation
    DVectorSlice::const_iterator
      d_itr = diag.begin(),
      x_itr = vs_rhs2.begin();
    DVectorSlice::iterator
      y_itr = vs_lhs->begin(),
      y_end = vs_lhs->end();
    for( ; y_itr != y_end; ++y_itr, ++d_itr, ++x_itr ) {
#ifdef LINALGPACK_CHECK_RANGE
      TEST_FOR_EXCEPT( !(  d_itr < diag.end()  ) );
      TEST_FOR_EXCEPT( !(  x_itr < vs_rhs2.end()  ) );
      TEST_FOR_EXCEPT( !(  y_itr < vs_lhs->end()  ) );
#endif
      *y_itr = beta * (*y_itr) + alpha * (*d_itr) * (*x_itr);
    }
  }
}

void MatrixSymDiagStd::Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  const DVectorSlice diag = this->diag();
  size_type n = diag.size();

  DenseLinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
    , n, n, trans_rhs1, sv_rhs2.size() );
  //
  // y = b*y + a * op(A) * x
  //
  DenseLinAlgPack::Vt_S(vs_lhs,beta); // y = b * y
  //
  // A is symmetric and diagonal A = diag(diag) so:
  //
  // y(j) += a * diag(j) * x(j), for j = 1...n
  //
  // x is sparse so take account of this.

  for(   SpVectorSlice::const_iterator x_itr = sv_rhs2.begin()
     ; x_itr != sv_rhs2.end()
     ; ++x_itr )
  {
    (*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
      += alpha * diag(x_itr->indice() + sv_rhs2.offset()) * x_itr->value();
      // Note: The indice x(i) invocations are ranged check
      // if this is compiled into the code.
  }
}

// Overridden from MatrixWithOpFactorized

void MatrixSymDiagStd::V_InvMtV(
  DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2) const
{
  const DVectorSlice diag = this->diag();
  size_type n = diag.size();

  // y = inv(op(A)) * x
  //
  // A is symmetric and diagonal (A = diag(diag)) so:
  //
  // y(j) = x(j) / diag(j), for j = 1...n

  DenseLinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
    , n, n, trans_rhs1, vs_rhs2.size() );
  
  if( vs_rhs2.stride() == 1 && vs_lhs->stride() == 1 ) {
    // Optimized implementation
    const value_type
      *d_itr      = diag.raw_ptr(),
      *x_itr      = vs_rhs2.raw_ptr();
    value_type
      *y_itr      = vs_lhs->raw_ptr(),
      *y_end      = y_itr + vs_lhs->size();
    while( y_itr != y_end )
      *y_itr++ = (*x_itr++) / (*d_itr++);
  }
  else {
    // Generic implementation
    DVectorSlice::const_iterator
      d_itr = diag.begin(),
      x_itr = vs_rhs2.begin();
    DVectorSlice::iterator
      y_itr = vs_lhs->begin(),
      y_end = vs_lhs->end();
    for( ; y_itr != y_end; ++y_itr, ++d_itr, ++x_itr ) {
      TEST_FOR_EXCEPT( !(  d_itr < diag.end()  ) );
      TEST_FOR_EXCEPT( !(  x_itr < vs_rhs2.end()  ) );
      TEST_FOR_EXCEPT( !(  y_itr < vs_lhs->end()  ) );
      *y_itr = (*x_itr)/(*d_itr);
    }
  }
}

void MatrixSymDiagStd::V_InvMtV(
  DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2) const
{
  const DVectorSlice diag = this->diag();
  size_type n = diag.size();

  // y = inv(op(A)) * x
  //
  // A is symmetric and diagonal A = diag(diag) so:
  //
  // y(j) = x(j) / diag(j), for j = 1...n
  //
  // x is sparse so take account of this.
  
  DenseLinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
    , n, n, trans_rhs1, sv_rhs2.size() );

  for(   SpVectorSlice::const_iterator x_itr = sv_rhs2.begin()
     ; x_itr != sv_rhs2.end()
     ; ++x_itr )
  {
    (*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
      = x_itr->value() / diag(x_itr->indice() + sv_rhs2.offset());
      // Note: The indice x(i) invocations are ranged check
      // if this is compiled into the code.
  }
}

} // end namespace AbstractLinAlgPack

#endif // 0
