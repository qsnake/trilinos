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
//

#ifndef LIN_ALG_OP_PACK_HACK_H
#define LIN_ALG_OP_PACK_HACK_H

#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp"

namespace LinAlgOpPack {

using DenseLinAlgPack::DVector;
using DenseLinAlgPack::DVectorSlice;
using DenseLinAlgPack::DMatrixSlice;
using AbstractLinAlgPack::SpVectorSlice;
using AbstractLinAlgPack::GenPermMatrixSlice;
using AbstractLinAlgPack::MatrixOp;
using AbstractLinAlgPack::MatrixNonsing;
using AbstractLinAlgPack::MatrixOpNonsing;

/** \brief gms_lhs = op(M_rhs). */
void assign(DMatrixSlice* gms_lhs, const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs);

/** \brief <tt>m_lhs += alpha * op(mwo_rhs1)</tt>.
 */
void Mp_StM(
  DMatrixSlice* vs_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  );

/** \brief <tt>vs_lhs = alpha * op(mwo_rhs1) * vs_rhs2 + beta * vs_lhs</tt>.
 */
void Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, const MatrixOp& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2
  ,value_type beta = 1.0 );

/** \brief <tt>vs_lhs = op(mwo_rhs1) * vs_rhs2 + beta * vs_lhs</tt>.
 */
void Vp_MtV(
  DVectorSlice* vs_lhs, const MatrixOp& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2
  ,value_type beta = 1.0 );

/** \brief <tt>vs_lhs = alpha * op(mwo_rhs1) * vs_rhs2 + beta * sv_lhs</tt>.
 */
void Vp_StMtV(
  DVectorSlice* vs_lhs, value_type alpha, const MatrixOp& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
  ,value_type beta = 1.0 );

/** \brief <tt>vs_lhs = op(mwo_rhs1) * vs_rhs2</tt>.
 */
void V_MtV(
  DVectorSlice* vs_lhs, const MatrixOp& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

/** \brief <tt>vs_lhs = inv(op(mwo_rhs1)) * vs_rhs2</tt>.
 */
void V_InvMtV(
  DVectorSlice* vs_lhs, const MatrixOpNonsing& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2 );

/** \brief <tt>v_lhs = inv(op(mwo_rhs1)) * vs_rhs2</tt>.
 */
void V_InvMtV(
  DVector* v_lhs, const MatrixOpNonsing& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2 );

/** \brief <tt>vs_lhs = inv(op(mwo_rhs1)) * sv_rhs2</tt>.
 */
void V_InvMtV(
  DVectorSlice* vs_lhs, const MatrixOpNonsing& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

/** \brief <tt>v_lhs = inv(op(mwo_rhs1)) * sv_rhs2</tt>.
 */
void V_InvMtV(
  DVector* v_lhs, const MatrixOpNonsing& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );

/** \brief <tt>vs_lhs = alpha * op(gpms_rhs1) * op(mwo_rhs2) * vs_rhs3 + beta * vs_lhs</tt>.
 */
void Vp_StPtMtV(
  DVectorSlice* vs_lhs, value_type alpha
  ,const GenPermMatrixSlice& gpms_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ,const DVectorSlice& vs_rhs3, value_type beta = 1.0 );

/** \brief <tt>vs_lhs = alpha * op(gpms_rhs1) * op(mwo_rhs2) * sv_rhs3 + beta * vs_lhs</tt>.
 */
void Vp_StPtMtV(
  DVectorSlice* vs_lhs, value_type alpha
  ,const GenPermMatrixSlice& gpms_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ,const SpVectorSlice& sv_rhs3, value_type beta = 1.0 );

} // end namespace LinAlgOpPack

// //////////////////////
// Inline functions

inline
void LinAlgOpPack::Vp_MtV(
  DVectorSlice* vs_lhs, const MatrixOp& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2
  ,value_type beta )
{
  Vp_StMtV(vs_lhs,1.0,mwo_rhs1,trans_rhs1,vs_rhs2,beta);
}

inline
void LinAlgOpPack::V_MtV(
  DVectorSlice* vs_lhs, const MatrixOp& mwo_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 )
{
  Vp_StMtV(vs_lhs,1.0,mwo_rhs1,trans_rhs1,sv_rhs2);
}




#endif // LIN_ALG_OP_PACK_HACK_H
