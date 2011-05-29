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

#ifndef MATRIX_HESSIAN_RELAXED_H
#define MATRIX_HESSIAN_RELAXED_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymOp.hpp"

namespace ConstrainedOptPack {

/** \brief Represents a symmetric Hessian matrix with a relaxation variable
  * added.
  *
  * This class is used to represent the matrix:
  \begin{verbatim}
       [ H       ]
  G =  [    bigM ]
  \end{verbatim}
  *
  */
class MatrixHessianRelaxed : public MatrixSymOp {
public:

  /// Construct to uninitialized
  MatrixHessianRelaxed();

  /** \brief Initialize.
    *
    * ToDo: Finish documentation!
    *
    */
  void initialize(
      const MatrixSymOp	&H
    , value_type			bigM
    );

  // ///////////////////////////////
  // Overridden from Matrix

  /** \brief . */
  size_type rows() const;

  // //////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const DVectorSlice& vs_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const SpVectorSlice& sv_rhs3, value_type beta) const;
  /** \brief . */
  value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const ;

private:
  size_type				n_;	// size of H
  const MatrixSymOp	*H_;
  value_type				bigM_;

};	// end class MatrixHessianRelaxed

}	// end namespace ConstrainedOptPack

#endif 	// MATRIX_HESSIAN_RELAXED_H
