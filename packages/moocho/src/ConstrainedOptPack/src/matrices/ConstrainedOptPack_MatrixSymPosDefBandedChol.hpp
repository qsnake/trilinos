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

#ifndef MATRIX_SYM_POS_DEF_BANDED_CHOL_H
#define MATRIX_SYM_POS_DEF_BANDED_CHOL_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/MatrixSymWithOpFactorized.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Miref_count_ptr.h"
#include "MiReleaseResource.h"

namespace ConstrainedOptPack {
/** \brief Matrix subclass for banded symmetric positive definite matrices and their
 * Cholesky factors.
 * 
 * This class is designed to support the LAPACK routines for banded symmetric positive
 * definite matrices.  The banded matrix and/or its cholesky factor are stored in
 * simple flat rectangular matrices compatible with the LAPACK routines.
 *
 * For example, for #n = 8, kd = 3# the original matrix #M# (if set) is stored in the
 * following format #MB#:
 \begin{verbatim}

         M                                   MB

 [ x x x x         ]                  [ x x x x x x x x ]
 [ x x x x x       ]  lower triangle  [ x x x x x x x o ]
 [ x x x x x x     ]         =>       [ x x x x x x o o ]
 [ x x x x x x x   ]                  [ x x x x x o o o ]
 [   x x x x x x x ]
 [     x x x x x x ]                  [ o o o x x x x x ]
 [       x x x x x ]  upper triangle  [ o o x x x x x x ]
 [         x x x x ]         =>       [ o x x x x x x x ]
                                      [ x x x x x x x x ]
 \end{verbatim}
 * The Cholesky factor #U# is sorted in a similar format #UB#.  Technically, the matrix
 * is #M = scale * U'*U# so that #M# may be negative definite as well.
 */
class MatrixSymPosDefBandedChol	: public MatrixSymWithOpFactorized
{
public:
  
  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>  release_resource_ptr_t;

  // //////////////
    // Constructors

  /** \brief Construct and Initialize.
   *
   * This constructor just calls #this->initialize(...)#.
   */
  MatrixSymPosDefBandedChol(
    size_type                         n                       = 0
    ,size_type                        kd                      = 0
    ,DMatrixSlice                   *MB                     = NULL
    ,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
    ,BLAS_Cpp::Uplo                   MB_uplo                 = BLAS_Cpp::lower
    ,DMatrixSlice                   *UB                     = NULL
    ,const release_resource_ptr_t&    UB_release_resource_ptr = NULL
    ,BLAS_Cpp::Uplo                   UB_uplo                 = BLAS_Cpp::lower
    ,bool                             update_factor           = false
    ,value_type                       scale                   = 1.0
    );

  // ///////////////////////////
    // Access representation

  /** \brief Initialize
   *
   * If called with all of the default arguments then #this# will become uninitialized.
   *
   * ToDo: Finish pre and post conditions!
   *
   * @param  n        [in] Determines the size of the banded matrix (n x n).
   *                  If #n == 0# then all of the following arguments should be left at
   *                  their defaults and #this# will become uninitialized.
   * @param  kd       [in] Determines the band width of the matrix as defined by xPBTRF(...).
   * @param  MB       [in/state] If #MB != NULL# then this matrix (size (kd+1) x n) is used to store
   *                  the original banded matrix #MB# in the format of xPBTRF(...).  This matrix must
   *                  be initialized on input.
   * @param  MB_release_resource_ptr
   *                  [in] Only significant if #MB != NULL#.  Points to a resource to
   *                  be released when #MB# is no longer needed.
   * @param  MB_uplo  [in] Determines if #MB# is stores the upper or lower triangular elements.
   * @param  UB       [in/state] If #UB != NULL# then this matrix (size (kd+1) x n) is used to store
   *                  the Cholesky factor of the banded matrix in the format of xPBTRF(...).
   *                  This matrix may or may not be initialized on input.
   *                  If #update_factor == false# this this matrix must  already be initialized.
   *                  If #update_factor == true# then this matrix will be computed.
2	 *                  If #UB == NULL# then storage for the Cholesky factor will be computed
   *                  on the fly and will be factored.
   * @param  UB_release_resource_ptr
   *                  [in] Only significant if #UB != NULL#.  Points to a resource to
   *                  be released when #UB# is no longer needed.
   * @param  UB_uplo  [in] Determines if #UB# is stores the upper or lower triangular elements.
   * @param  update_factor
   *                  [in] If true then the factor will be computed within this function call.
   * @param  scale    [in] Only significant if #MB != NULL# or #UB != NULL# (see intro).
   */
  void initialize(
    size_type                         n                       = 0
    ,size_type                        kd                      = 0
    ,DMatrixSlice                   *MB                     = NULL
    ,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
    ,BLAS_Cpp::Uplo                   MB_uplo                 = BLAS_Cpp::lower
    ,DMatrixSlice                   *UB                     = NULL
    ,const release_resource_ptr_t&    UB_release_resource_ptr = NULL
    ,BLAS_Cpp::Uplo                   UB_uplo                 = BLAS_Cpp::lower
    ,bool                             update_factor           = false
    ,value_type                       scale                   = 1.0
    );

  /** \brief . */
  size_type kd() const;
  /** \brief Get view of MB.
   */
  DMatrixSlice& MB();
  /** \brief . */
  const DMatrixSlice& MB() const;
  /** \brief . */
  BLAS_Cpp::Uplo MB_uplo() const;
  /** \brief Get view of UB.
   */
  DMatrixSlice& UB();
  /** \brief . */
  const DMatrixSlice& UB() const;
  /** \brief . */
  BLAS_Cpp::Uplo UB_uplo() const;

  // /////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type nz() const;
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
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

  // //////////////////////////////////
  // Overridden from MatrixFactorized

  /// With throw exception if factorization is not allowed.
  void V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2) const;

private:
  
  // /////////////////////////////
  // Private data members

  size_type                       n_;
  size_type                       kd_;
  DMatrixSlice                  MB_;
  release_resource_ptr_t          MB_release_resource_ptr_;
  BLAS_Cpp::Uplo                  MB_uplo_;
  mutable DMatrixSlice          UB_;
  mutable release_resource_ptr_t  UB_release_resource_ptr_;
  mutable BLAS_Cpp::Uplo          UB_uplo_;
  mutable bool                    factor_updated_;
  value_type                      scale_;

  // /////////////////////////////
  // Private member functions

  void assert_initialized() const;
  void update_factorization() const;

}; // end class MatrixSymPosDefBandedChol

// ///////////////////////////////////////////////////////
// Inline members for MatrixSymPosDefBandedChol

inline
size_type MatrixSymPosDefBandedChol::kd() const
{
  return kd_;
}

inline
DMatrixSlice& MatrixSymPosDefBandedChol::MB()
{
  return MB_;
}

inline
const DMatrixSlice& MatrixSymPosDefBandedChol::MB() const
{
  return MB_;
}

inline
BLAS_Cpp::Uplo MatrixSymPosDefBandedChol::MB_uplo() const
{
  return MB_uplo_;
}

inline
DMatrixSlice& MatrixSymPosDefBandedChol::UB()
{
  return UB_;
}

inline
const DMatrixSlice& MatrixSymPosDefBandedChol::UB() const
{
  return UB_;
}

inline
BLAS_Cpp::Uplo MatrixSymPosDefBandedChol::UB_uplo() const
{
  return UB_uplo_;
}

} // end namespace ConstrainedOptPack

#endif // MATRIX_SYM_POS_DEF_BANDED_CHOL_H
