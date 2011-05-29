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

#include <iostream> // Debuggin only

#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack {

MatrixSymDiagStd::MatrixSymDiagStd(
  const VectorSpace::vec_mut_ptr_t& diag
  ,bool                             unique
  )
{
  this->initialize(diag,unique);
//	std::cerr << "MatrixSymDiagStd::rows() = " << this->rows() << std::endl; // Debugging
//	std::cerr << "MatrixSymDiagStd::nz() = "   << this->nz()   << std::endl; // Debugging
//	std::cerr << "MatrixSymDiagStd::cols() = " << this->cols() << std::endl; // Debugging
//	std::cerr << "MatrixSymDiagStd::nz() = "   << this->nz()   << std::endl; // Debugging
}

void MatrixSymDiagStd::initialize(
  const VectorSpace::vec_mut_ptr_t& diag
  ,bool                             unique
  )
{
  diag_   = diag;   // lazy copy!
  unique_ = unique;
}

VectorMutable& MatrixSymDiagStd::diag()
{
  copy_unique();
  VectorMutable *diag = diag_.get();
  TEST_FOR_EXCEPTION(
    !diag, std::logic_error
    ,"MatrixSymDiagStd::diag(): Error, the diagonal vector has not been set! " );
  return *diag;;
}

const VectorSpace::vec_mut_ptr_t&
MatrixSymDiagStd::diag_ptr() const
{
  return diag_;
}

// Overridden from MatrixBase

size_type MatrixSymDiagStd::rows() const
{
  return diag_.get() ? diag_->dim() : 0;
}

size_type MatrixSymDiagStd::nz() const
{
  return diag_.get() ? diag_->nz() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixSymDiagStd::space_cols() const {
  return diag_->space();
}

const VectorSpace& MatrixSymDiagStd::space_rows() const {
  return diag_->space();
}

MatrixOp&
MatrixSymDiagStd::operator=(const MatrixOp& M)
{
  const MatrixSymDiagStd
    *p_M = dynamic_cast<const MatrixSymDiagStd*>(&M);

  TEST_FOR_EXCEPTION(
    p_M == NULL, std::logic_error
    ,"MatrixSymDiagStd::operator=(M): Error, the matrix M with concrete type "
    "\'" << typeName(M) << "\' does not support the MatrixSymDiagStd type! " );

  if( p_M == this ) return *this; // Assignment to self

  diag_    = p_M->diag_;  // lazy copy!
  unique_  = p_M->unique_;

  return *this;
}

bool MatrixSymDiagStd::Mp_StM(
  MatrixOp* M_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const
{
  // ToDo: validate the vector spaces for the matrices!
  MultiVectorMutable
    *M_mv_lhs = dynamic_cast<MultiVectorMutable*>(M_lhs);
  if(!M_mv_lhs)
    return false;
  VectorSpace::vec_mut_ptr_t
    M_diag = M_mv_lhs->diag(0);
  if(!M_diag.get())
    return false; // Access to the diagonal is not supported!
  Vp_StV( M_diag.get(), alpha, *diag_ );
  return true;
}

void MatrixSymDiagStd::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const Vector& v_rhs2, value_type beta) const
{
  // ToDo: Validate input!
  if(beta == 0.0)
    *v_lhs = 0.0;
  else if(beta != 1.0)
    Vt_S( v_lhs, beta );
  ele_wise_prod( alpha, v_rhs2, *diag_, v_lhs );
}

void MatrixSymDiagStd::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  MatrixOp::Vp_StMtV(v_lhs,alpha,trans_rhs1,sv_rhs2,beta); // ToDo: Implement specialized!
}

// Overridden from MatrixNonsing

void MatrixSymDiagStd::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  , const Vector& v_rhs2) const
{
  ele_wise_divide( 1.0, v_rhs2, *diag_, v_lhs );
}

void MatrixSymDiagStd::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2) const
{
  MatrixNonsing::V_InvMtV(v_lhs,trans_rhs1,sv_rhs2 ); // ToDo: Implement specialized!
}

bool MatrixSymDiagStd::syrk(
  BLAS_Cpp::Transp   A_trans
  ,value_type        a
  ,value_type        b
  ,MatrixSymOp   *B
  ) const
{
  MatrixSymDiagStd    *B_sd = dynamic_cast<MatrixSymDiagStd*>(B);
  if(!B_sd) return false;
  VectorMutable     &B_diag = B_sd->diag();
  const Vector      &A_diag = this->diag();
  // B = b*B + a*A*A
  Vt_S( &B_diag, b );
  ele_wise_prod( 1.0, A_diag, A_diag, &B_diag );   // B.diag(i) += a * (A.diag)(i) * (A.diag)(i)
  return true;
}

// Overridden from MatrixSymInitDiag

void MatrixSymDiagStd::init_identity( const VectorSpace& space_diag, value_type alpha )
{
  diag_ = space_diag.create_member();
  if( diag_->dim() )
    *diag_ = alpha;
}

void MatrixSymDiagStd::init_diagonal( const Vector& diag )
{
  diag_ = diag.space().create_member();
  *diag_ = diag;
}

// Overridden from MatrixSymDiag

const Vector& MatrixSymDiagStd::diag() const
{
  return const_cast<MatrixSymDiagStd*>(this)->diag();
}

// private

void MatrixSymDiagStd::copy_unique()
{
  if( diag_.get() && diag_.count() > 1 && unique_ )
    diag_ = diag_->clone();
}

} // end namespace AbstractLinAlgPack
