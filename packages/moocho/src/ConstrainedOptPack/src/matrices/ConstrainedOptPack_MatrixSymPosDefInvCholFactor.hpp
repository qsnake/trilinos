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

#ifndef MATRIX_SYM_POS_DEF_INV_CHOL_FACTOR_H
#define MATRIX_SYM_POS_DEF_INV_CHOL_FACTOR_H

#include "SymInvCholMatrixClass.hpp"
#include "AbstractLinAlgPack_MatrixSymSecant.hpp"
#include "AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp"
#include "AbstractLinAlgPack_MatrixWithOpConcreteEncap.hpp"
#include "AbstractLinAlgPack/src/MatrixSymWithOpFactorized.hpp"
#include "SerializationPack_Serializable.hpp"

namespace ConstrainedOptPack {

/** \brief Implementation of MatrixOp abstract interface for SymInvCholMatrix
  */
class MatrixSymPosDefInvCholFactor
  : public virtual MatrixWithOpConcreteEncap<SymInvCholMatrix>
  , public virtual MatrixSymWithOpFactorized
  , public MatrixSymSecant
  , public MatrixExtractInvCholFactor
  , public virtual SerializationPack::Serializable
{
public:

  /** \brief . */
  MatrixSymPosDefInvCholFactor()
  {}

  /** \brief . */
  MatrixSymPosDefInvCholFactor(const SymInvCholMatrix& m)
    : MatrixWithOpConcreteEncap<SymInvCholMatrix>(m)
  {}

  /** @name Overridden from Matrix */
  //@{

  /// 
  size_type cols() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  MatrixOp& operator=(const MatrixOp& m);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  value_type transVtMtV(const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const DVectorSlice& vs_rhs3) const;
  /** \brief . */
  value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const;

  //@}

  /** @name Overridden from MatrixFactorized */
  //@{

  /** \brief . */
  void V_InvMtV(DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2) const;
  /** \brief . */
  void V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2) const;
  /** \brief . */
  void V_InvMtV(DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2) const;
  /** \brief . */
  void V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2) const;
  /** \brief . */
  value_type transVtInvMtV(const DVectorSlice& vs_rhs1
    , BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3) const;
  /** \brief . */
  value_type transVtInvMtV(const SpVectorSlice& sv_rhs1
    , BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;

  //@}

  /** @name Overridden from MatrixSymFactorized */
  //@{

  /** \brief . */
  void M_StMtInvMtM( DMatrixSliceSym* sym_gms_lhs, value_type alpha
    , const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans, EMatrixDummyArg
    ) const;

  //@}

  /** @name Overridden from MatrixSymSecant */
  //@{

  /** \brief . */
  void init_identity( size_type n, value_type alpha );
  /** \brief . */
  void init_diagonal( const DVectorSlice& diag );
  /** \brief . */
  void secant_update(DVectorSlice* s, DVectorSlice* y, DVectorSlice* Bs);

  //@}

  /** @name Overridden from MatrixExtractInvCholFactor */
  //@{

  /** \brief . */
  void extract_inv_chol( DMatrixSliceTriEle* InvChol ) const;

  //@}

  /** @name Overridden from Serializable */
  //@{

  /** \brief . */
  void serialize( std::ostream &out ) const;
  /** \brief . */
  void unserialize( std::istream &in );

  //@}

};	// end class MatrixSymPosDefInvCholFactor

}	// end namespace ConstrainedOptPack 

#endif	// MATRIX_SYM_POS_DEF_INV_CHOL_FACTOR_H
