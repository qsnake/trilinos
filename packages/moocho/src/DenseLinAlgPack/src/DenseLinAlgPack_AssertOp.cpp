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

#include <stdexcept>
#include <string>

#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef LINALGPACK_CHECK_RHS_SIZES

void DenseLinAlgPack::Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size)
{
  TEST_FOR_EXCEPTION(
    v_lhs_size != v_rhs_size, std::length_error
    ,"Vp_V_assert_sizes(...) : The sizes of v_lhs = " << v_lhs_size << " and v_rhs = " << v_rhs_size
    << " in the operation v_lhs += op v_rhs do not match");
}

void DenseLinAlgPack::VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size)
{
  TEST_FOR_EXCEPTION(
    v_rhs1_size != v_rhs2_size, std::length_error
    ,"VopV_assert_sizes(...) : The sizes of v_rhs1 and v_rhs2 "
    "in the operation v_rhs1 op v_rhs2 do not match");
}

void DenseLinAlgPack::Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs)
{
  if(		rows(m_lhs_rows,m_lhs_cols,trans_lhs) != rows(m_rhs_rows,m_rhs_cols,trans_rhs)
    ||	cols(m_lhs_rows,m_lhs_cols,trans_lhs) != cols(m_rhs_rows,m_rhs_cols,trans_rhs) )
  {
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_M_assert_sizes(...) : The sizes of m_lhs and m_rhs "
      "in the operation op(m_lhs) += op op(m_rhs) do not match");
  }
}

void DenseLinAlgPack::MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
  if(		rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2)
    ||	cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != cols(m_rhs2_rows,m_rhs2_cols,trans_rhs2) )
  {
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_M_assert_sizes(...) : The sizes of m_rhs1 and m_rhs2 "
      "in the operation op(m_rhs1) op op(m_rhs2) do not match");
  }
}

void DenseLinAlgPack::MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
  , BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_rhs2_size)
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"MtV_assert_sizes(...) : The number of columns in "
      "m_rhs1 and the size of v_rhs2 in the operation v_lhs += op(m_rhs1) * v_rhs2 "
      "do not match");
}

void DenseLinAlgPack::Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
  , size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_rhs2_size)
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Vp_MtV_assert_sizes(...) : The number of columns in"
      " m_rhs1 and the size of v_rhs2 in the operation v_lhs += op(m_rhs1) * v_rhs2"
      " do not match");
  if(rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_lhs_size)
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Vp_MtV_assert_sizes(...) : The number of rows in"
      " m_rhs1 and the size of v_lhs in the operation v_lhs += op(m_rhs1) * v_rhs2"
      " do not match");
}

void DenseLinAlgPack::MtM_assert_sizes(
    size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"MtM_assert_sizes(...) : The number of columns in"
      " m_rhs1 and the number of rows in m_rhs2 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
}

void DenseLinAlgPack::Mp_MtM_assert_sizes(
    size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_MtM_assert_sizes(...) : The number of columns in"
      " m_rhs1 and the number of rows in m_rhs2 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
  if(rows(m_lhs_rows,m_lhs_cols,trans_lhs) != rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1))
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_MtM_assert_sizes(...) : The number of rows in"
      " m_lhs and the number of rows in m_rhs1 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
  if(cols(m_lhs_rows,m_lhs_cols,trans_lhs) != cols(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
    TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_MtM_assert_sizes(...) : The number of columns in"
      " m_lhs and the number of columns in m_rhs1 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
}

#endif // LINALGPACK_CHECK_RHS_SIZES
