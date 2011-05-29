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

#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

// Level 1 BLAS for Matrices

// M_lhs = op(M_rhs).

void LinAlgOpPack::assign(MatrixOp* M_lhs, const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
  Mp_M_assert_compatibility( M_lhs, BLAS_Cpp::no_trans, M_rhs, trans_rhs );
  M_lhs->zero_out();
  Mp_StM(M_lhs,1.0,M_rhs,trans_rhs);
}

// M_lhs = alpha * op(M_rhs).

void LinAlgOpPack::M_StM(MatrixOp* M_lhs, value_type alpha, const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
  Mp_M_assert_compatibility( M_lhs, BLAS_Cpp::no_trans, M_rhs, trans_rhs );
  M_lhs->zero_out();
  Mp_StM(M_lhs,alpha,M_rhs,trans_rhs);
}

// M_lhs = op(M_rhs1) + op(M_rhs2).

void LinAlgOpPack::M_MpM(MatrixOp* M_lhs, const MatrixOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
  MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  M_lhs->zero_out();
  Mp_M(M_lhs,M_rhs1,trans_rhs1);
  Mp_M(M_lhs,M_rhs2,trans_rhs2);
}

// M_lhs = op(M_rhs1) - op(M_rhs2).

void LinAlgOpPack::M_MmM(MatrixOp* M_lhs, const MatrixOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
  MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  M_lhs->zero_out();
  Mp_M(M_lhs,M_rhs1,trans_rhs1);
  Mp_StM(M_lhs,-1.0,M_rhs2,trans_rhs2);
}

// M_lhs = alpha * op(M_rhs1) + op(m_rhs2).

void LinAlgOpPack::M_StMpM(MatrixOp* M_lhs, value_type alpha, const MatrixOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
  MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  assign(M_lhs,M_rhs2,trans_rhs2);
  Mp_StM(M_lhs,alpha,M_rhs1,trans_rhs1);
}

// Level 3 BLAS

// M_lhs = alpha * op(M_rhs1) * op(M_rhs2).

void LinAlgOpPack::M_StMtM(MatrixOp* M_lhs, value_type alpha, const MatrixOp& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_MtM_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  Mp_StMtM(M_lhs,alpha,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,0.0);
}

// M_lhs = op(M_rhs1) * op(M_rhs2).

void LinAlgOpPack::M_MtM(MatrixOp* M_lhs, const MatrixOp& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_MtM_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,0.0);
}
