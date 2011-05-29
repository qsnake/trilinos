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

#include "AbstractLinAlgPack_MatrixZero.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack {

// Constructors/initializers

MatrixZero::MatrixZero(
  const VectorSpace::space_ptr_t&    space_cols
  ,const VectorSpace::space_ptr_t&   space_rows
  )
{
  this->initialize(space_cols,space_rows);
}

void MatrixZero::initialize(
  const VectorSpace::space_ptr_t&    space_cols
  ,const VectorSpace::space_ptr_t&   space_rows
  )
{
  TEST_FOR_EXCEPTION(
    (space_cols.get() == NULL && space_rows.get() != NULL)
    || (space_cols.get() != NULL && space_rows.get() == NULL)
    , std::invalid_argument
    ,"MatrixZero::initialize(...) : Error, the space_cols.get() and "
    "space_rows.get() must both be != NULL or == NULL" );
  space_cols_ = space_cols;
  space_rows_ = space_rows;
}

// Overridden from MatrixBase

size_type MatrixZero::rows() const
{
  return space_cols_.get() ? space_cols_->dim() : 0;
}

size_type MatrixZero::cols() const
{
  return space_rows_.get() ? space_rows_->dim() : 0;
}

size_type MatrixZero::nz() const
{
  return 0;
}

// Overridden form MatrixOp

const VectorSpace& MatrixZero::space_cols() const
{
  assert_initialized();
  return *space_cols_;
}

const VectorSpace& MatrixZero::space_rows() const
{
  assert_initialized();
  return *space_rows_;
}

void MatrixZero::zero_out()
{
  assert_initialized();
  // Automatically satisfied!
}

void MatrixZero::Mt_S( value_type alpha )
{
  assert_initialized();
  // Automatically satisfied!
}

MatrixOp& MatrixZero::operator=(const MatrixOp& M)
{
  assert_initialized();
  TEST_FOR_EXCEPT(true); // ToDo: Implement!
  return *this;
}

std::ostream& MatrixZero::output(std::ostream& out) const
{
  assert_initialized();
  return out << "Zero matrix of dimension " << rows() << " x " << cols() << std::endl;
}

// Level-1 BLAS

bool MatrixZero::Mp_StM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  assert_initialized();
  return true; // Nothing to do!
}

bool MatrixZero::Mp_StMtP(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  assert_initialized();
  return true; // Nothing to do!
}

bool MatrixZero::Mp_StPtM(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  assert_initialized();
  return true; // Nothing to do!
}

bool MatrixZero::Mp_StPtMtP(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  assert_initialized();
  return true; // Nothing to do!
}

// Level-2 BLAS

void MatrixZero::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans_in
  , const Vector& x, value_type b
  ) const
{
  assert_initialized();
  Vt_S(y,b);
}

void MatrixZero::Vp_StMtV(
  VectorMutable* y, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& x, value_type b) const
{
  assert_initialized();
  Vt_S(y,b);
}

void MatrixZero::Vp_StPtMtV(
  VectorMutable* y, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const Vector& x, value_type b) const
{
  assert_initialized();
  Vt_S(y,b);
}

void MatrixZero::Vp_StPtMtV(
  VectorMutable* y, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const SpVectorSlice& x, value_type b) const
{
  assert_initialized();
  Vt_S(y,b);
}

value_type MatrixZero::transVtMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const Vector& v_rhs3) const
{
  assert_initialized();
  return 0.0; // Nothing to do!
}

value_type MatrixZero::transVtMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const SpVectorSlice& sv_rhs3) const
{
  assert_initialized();
  return 0.0; // Nothing to do!
}

void MatrixZero::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  assert_initialized();
  sym_lhs->Mt_S(beta);
}

// Level-3 BLAS

bool MatrixZero::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  assert_initialized();
  m_lhs->Mt_S(beta);
  return true;
}

bool MatrixZero::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  , const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  , BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
  assert_initialized();
  m_lhs->Mt_S(beta);
  return true;
}

bool MatrixZero::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  assert_initialized();
  sym_lhs->Mt_S(beta);
  return true;
}

// private

void MatrixZero::assert_initialized() const {
  TEST_FOR_EXCEPTION(
    space_cols_.get() == NULL, std::logic_error
    ,"Error, the MatrixZero object has not been initialized!" );
}

} // end namespace AbstractLinAlgPack
