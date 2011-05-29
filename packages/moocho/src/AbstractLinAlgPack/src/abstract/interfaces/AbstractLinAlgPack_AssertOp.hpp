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

#ifndef ABSTRACT_LIN_ALG_PACK_ASSERT_OP_H
#define ABSTRACT_LIN_ALG_PACK_ASSERT_OP_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \defgroup AbstractLinAlgPackAssertOp_funcs Assertion functions for linear algebra operations.
  *
  * These functions check the compatibility of the vector spaces for many different types
  * of linear algebra operations and throw <tt>VectorSpace::IncompatibleVectorSpaces</tt> 
  * expressions if the vector spaces do not match.  These functions only perform there checks
  * if <tt>ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY</tt> is defined.  These functions will also throw
  * <tt>std::invalid_argument</tt> if a lhs argument is <tt>NULL</tt>.
  */
//@{

/** @name Level 1 BLAS
  */
//@{
/// v_lhs += op v_rhs
void Vp_V_assert_compatibility(VectorMutable* v_lhs, const Vector& v_rhs);
/// v_lhs += op sv_rhs
void Vp_V_assert_compatibility(VectorMutable* v_lhs, const SpVectorSlice& sv_rhs);
/// v_rhs1 op v_rhs2
void VopV_assert_compatibility(const Vector& v_rhs1, const Vector& v_rhs2);
/// v_rhs1 op sv_rhs2
void VopV_assert_compatibility(const Vector& v_rhs1, const SpVectorSlice& sv_rhs2);
/// sv_rhs1 op v_rhs2
void VopV_assert_compatibility(const SpVectorSlice& sv_rhs1, const Vector& v_rhs2);
/// op(m_lhs) += op op(m_rhs)
void Mp_M_assert_compatibility(
  MatrixOp* m_lhs, BLAS_Cpp::Transp trans_lhs
  ,const MatrixOp& m_rhs, BLAS_Cpp::Transp trans_rhs );
/// op(m_rhs1) op op(m_rhs2)
void MopM_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 );
//@}

/**  @name Level 2 BLAS
  */
//@{
/// op(m_rhs1) * v_rhs2
void MtV_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2 );
/// op(m_rhs1) * sv_rhs2
void MtV_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );
/// v_lhs += op(m_rhs1) * v_rhs2
void Vp_MtV_assert_compatibility(
  VectorMutable* v_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2 );
/// v_lhs += op(m_rhs1) * sv_rhs2
void Vp_MtV_assert_compatibility(
  VectorMutable* v_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 );
//@}

/**  @name Level 3 BLAS
  */
//@{
/// op(m_lhs) += op(m_rhs1)
void MtM_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 );
/// op(m_lhs) += op(m_rhs1) * op(m_rhs2)
void Mp_MtM_assert_compatibility(
  MatrixOp* m_lhs, BLAS_Cpp::Transp trans_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 );
//@}

//@}

#ifndef ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY

// inline definitions that do nothing

inline
void Vp_V_assert_compatibility(VectorMutable* v_lhs, const Vector& v_rhs)
{} 
inline
void Vp_V_assert_compatibility(VectorMutable* v_lhs, const SpVectorSlice& sv_rhs)
{} 
inline
void VopV_assert_compatibility(const Vector& v_rhs1, const Vector&  v_rhs2)
{}
inline
void VopV_assert_compatibility(const Vector& v_rhs1, const SpVectorSlice& sv_rhs2)
{} 
inline
void VopV_assert_compatibility(const SpVectorSlice& sv_rhs1, const Vector& v_rhs2)
{} 
inline
void Mp_M_assert_compatibility(
  MatrixOp* m_lhs, BLAS_Cpp::Transp trans_lhs
  ,const MatrixOp& m_rhs, BLAS_Cpp::Transp trans_rhs )
{}
inline
void MopM_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{}
inline
void MtV_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2 )
{}
inline
void MtV_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 )
{} 
inline
void Vp_MtV_assert_compatibility(
  VectorMutable* v_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2 )
{}
inline
void Vp_MtV_assert_compatibility(
  VectorMutable* v_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 )
{}
inline
void MtM_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{}
inline
void Mp_MtM_assert_compatibility(
  MatrixOp* m_lhs, BLAS_Cpp::Transp trans_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{}

#endif // ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY

} // end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LIN_ALG_PACK_ASSERT_OP_H
