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

#include "AbstractLinAlgPack_MatrixOpNonsingAggr.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

MatrixOpNonsingAggr::MatrixOpNonsingAggr()
{} // Nothing to explicitly initialize

MatrixOpNonsingAggr::MatrixOpNonsingAggr(
  const mwo_ptr_t       &mwo
  ,BLAS_Cpp::Transp     mwo_trans
  ,const mns_ptr_t      &mns
  ,BLAS_Cpp::Transp     mns_trans
  )
{
  this->initialize(mwo,mwo_trans,mns,mns_trans);
}

void MatrixOpNonsingAggr::initialize(
  const mwo_ptr_t       &mwo
  ,BLAS_Cpp::Transp     mwo_trans
  ,const mns_ptr_t      &mns
  ,BLAS_Cpp::Transp     mns_trans
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    mwo.get() == NULL, std::invalid_argument
    ,"MatrixOpNonsingAggr::initialize(...): Error!" );
  TEST_FOR_EXCEPTION(
    mns.get() == NULL, std::invalid_argument
    ,"MatrixOpNonsingAggr::initialize(...): Error!" );
  const size_type
    mwo_rows = mwo->rows(),
    mwo_cols = mwo->cols(),
    mns_rows = mns->rows(),
    mns_cols = mns->cols();
  TEST_FOR_EXCEPTION(
    mwo_rows != mwo_cols, std::invalid_argument
    ,"MatrixOpNonsingAggr::initialize(...): Error!" );
  TEST_FOR_EXCEPTION(
    mns_rows != mns_cols, std::invalid_argument
    ,"MatrixOpNonsingAggr::initialize(...): Error!" );
  TEST_FOR_EXCEPTION(
    mwo_rows != mns_rows, std::invalid_argument
    ,"MatrixOpNonsingAggr::initialize(...): Error!" );
#endif
  mwo_       = mwo;
  mwo_trans_ = mwo_trans;
  mns_       = mns;
  mns_trans_ = mns_trans;
}

void MatrixOpNonsingAggr::set_uninitialized()
{
  namespace rcp = MemMngPack;
  mwo_       = Teuchos::null;
  mwo_trans_ = BLAS_Cpp::no_trans;
  mns_       = Teuchos::null;
  mns_trans_ = BLAS_Cpp::no_trans;
}

// Overridden from MatrixBase

size_type MatrixOpNonsingAggr::rows() const
{
  return mwo_.get() ? mwo_->rows() : 0; // square matrix!
}

size_type MatrixOpNonsingAggr::cols() const
{
  return mwo_.get() ? mwo_->rows() : 0; // square matrix!
}

size_type MatrixOpNonsingAggr::nz() const
{
  return mwo_.get() ? mwo_->nz() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixOpNonsingAggr::space_cols() const
{
  return mwo_trans_ == BLAS_Cpp::no_trans ? mwo_->space_cols() : mwo_->space_rows();
}

const VectorSpace& MatrixOpNonsingAggr::space_rows() const
{
  return mwo_trans_ == BLAS_Cpp::no_trans ? mwo_->space_rows() : mwo_->space_cols();
}

MatrixOp::mat_ptr_t
MatrixOpNonsingAggr::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  return MatrixOp::sub_view(row_rng,col_rng); // ToDo: Speicalize!
}

MatrixOp& MatrixOpNonsingAggr::operator=(const MatrixOp& M)
{
  using Teuchos::dyn_cast;
  const MatrixOpNonsingAggr
    Mp = dyn_cast<const MatrixOpNonsingAggr>(M);
  if( this == &Mp )
    return *this; // Assignment to self
  // Shallow copy is okay as long as client is careful!
  mwo_       = Mp.mwo_;
  mwo_trans_ = Mp.mwo_trans_;
  mns_       = Mp.mns_;
  mns_trans_ = Mp.mns_trans_;
  return *this;
}

std::ostream& MatrixOpNonsingAggr::output(std::ostream& out) const
{
  out << "Aggregate nonsingular matrix:\n";
  out << (mwo_trans_ == BLAS_Cpp::no_trans ? "mwo" : "mwo\'") << " =\n" << *mwo_;
  out << (mns_trans_ == BLAS_Cpp::no_trans ? "mns" : "mns\'") << " = ???\n";
  return out;
}

bool MatrixOpNonsingAggr::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  AbstractLinAlgPack::Mp_StM(mwo_lhs,alpha,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs));
  return true;
}

bool MatrixOpNonsingAggr::Mp_StMtP(
  MatrixOp* mwo_lhs, value_type alpha
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  AbstractLinAlgPack::Mp_StMtP(
    mwo_lhs,alpha,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans)
    ,P_rhs,P_rhs_trans);
  return true;
}

bool MatrixOpNonsingAggr::Mp_StPtM(
  MatrixOp* mwo_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  AbstractLinAlgPack::Mp_StPtM(
    mwo_lhs,alpha,P_rhs,P_rhs_trans
    ,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans));
  return true;
}

bool MatrixOpNonsingAggr::Mp_StPtMtP(
  MatrixOp* mwo_lhs, value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  AbstractLinAlgPack::Mp_StPtMtP(
    mwo_lhs,alpha,P_rhs1,P_rhs1_trans
    ,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans)
    ,P_rhs2,P_rhs2_trans);
  return true;
}

void MatrixOpNonsingAggr::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const Vector& x, value_type b) const
{
  AbstractLinAlgPack::Vp_StMtV(y,a,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),x,b);
}

void MatrixOpNonsingAggr::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b) const
{
  AbstractLinAlgPack::Vp_StMtV(y,a,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),x,b);
}

void MatrixOpNonsingAggr::Vp_StPtMtV(
  VectorMutable* vs_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const Vector& v_rhs3, value_type beta) const
{
  AbstractLinAlgPack::Vp_StPtMtV(
    vs_lhs,alpha,P_rhs1,P_rhs1_trans
    ,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_rhs2_trans)
    ,v_rhs3,beta);
}

void MatrixOpNonsingAggr::Vp_StPtMtV(
  VectorMutable* vs_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const SpVectorSlice& sv_rhs3, value_type beta) const
{
  AbstractLinAlgPack::Vp_StPtMtV(
    vs_lhs,alpha,P_rhs1,P_rhs1_trans
    ,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_rhs2_trans)
    ,sv_rhs3,beta);
}

value_type MatrixOpNonsingAggr::transVtMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const Vector& v_rhs3) const
{
  return AbstractLinAlgPack::transVtMtV(v_rhs1,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),v_rhs3);
}

value_type MatrixOpNonsingAggr::transVtMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  ,const SpVectorSlice& sv_rhs3
  ) const
{
  return AbstractLinAlgPack::transVtMtV(sv_rhs1,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),sv_rhs3);
}

void MatrixOpNonsingAggr::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  ,const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  ,const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  ,value_type beta, MatrixSymOp* symwo_lhs
  ) const
{
  AbstractLinAlgPack::syr2k(
    *mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans)
    ,alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs);
}

bool MatrixOpNonsingAggr::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
  ,BLAS_Cpp::Transp trans_rhs2, value_type beta
  ) const
{
  AbstractLinAlgPack::Mp_StMtM(
    mwo_lhs,alpha,*mwo_,trans_rhs1
    ,mwo_rhs2,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),beta);
  return true;
}

bool MatrixOpNonsingAggr::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2, value_type beta
  ) const
{
  AbstractLinAlgPack::Mp_StMtM(
    mwo_lhs,alpha,mwo_rhs1,trans_rhs1
    ,*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,trans_rhs2),beta);
  return true;
}

bool MatrixOpNonsingAggr::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  ,value_type beta, MatrixSymOp* sym_lhs
  ) const
{
  AbstractLinAlgPack::syrk(*mwo_,BLAS_Cpp::trans_trans(mwo_trans_,M_trans),alpha,beta,sym_lhs);
  return true;
}

// Overridden from MatrixNonsing */

void MatrixOpNonsingAggr::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& v_rhs2
  ) const
{
  AbstractLinAlgPack::V_InvMtV(v_lhs,*mns_,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1),v_rhs2);
}

void MatrixOpNonsingAggr::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  ,const SpVectorSlice& sv_rhs2
  ) const
{
  AbstractLinAlgPack::V_InvMtV(v_lhs,*mns_,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1),sv_rhs2);
}

value_type MatrixOpNonsingAggr::transVtInvMtV(
  const Vector& v_rhs1
  ,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3
  ) const
{
  return AbstractLinAlgPack::transVtInvMtV(v_rhs1,*mns_,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs2),v_rhs3);
}

value_type MatrixOpNonsingAggr::transVtInvMtV(
  const SpVectorSlice& sv_rhs1
  ,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3
  ) const
{
  return AbstractLinAlgPack::transVtInvMtV(sv_rhs1,*mns_,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs2),sv_rhs3);
}

void MatrixOpNonsingAggr::M_StInvMtM(
  MatrixOp* m_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ) const
{
  AbstractLinAlgPack::M_StInvMtM(m_lhs,alpha,*mns_,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1),mwo_rhs2,trans_rhs2);
}

void MatrixOpNonsingAggr::M_StMtInvM(
  MatrixOp* m_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  AbstractLinAlgPack::M_StMtInvM(m_lhs,alpha,mwo_rhs1,trans_rhs1,*mns_,BLAS_Cpp::trans_trans(mns_trans_,trans_rhs1));
}

} // end namespace AbstractLinAlgPack
