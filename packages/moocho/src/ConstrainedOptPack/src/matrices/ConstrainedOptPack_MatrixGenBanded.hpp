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

#ifndef MATRIX_GEN_BANDED_H
#define MATRIX_GEN_BANDED_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Miref_count_ptr.h"
#include "MiReleaseResource.h"

namespace ConstrainedOptPack {
/** \brief Matrix subclass for general (possibly singular) banded matrices.
 * 
 * The banded matrix is stored by column in a simple flat rectangular matrix.
 * For example, for #m = 10, n = 8, kl = 3, ku = 2# the matrix #M# is stored in the
 * following format #MB# (same as for the BLAS routine xGBMV(...)):
 \begin{verbatim}

         M                                   MB
 [ x x x               ]
 [ x x x x             ]         [ o o x x x x x x x x ] \ ku = 2
 [ x x x x x           ]         [ o x x x x x x x x o ] /
 [ x x x x x x         ]    =>   [ x x x x x x x x o o ]
 [   x x x x x x       ]         [ x x x x x x x o o o ] \
 [     x x x x x x     ]         [ x x x x x x o o o o ] | kl = 3
 [       x x x x x x   ]         [ x x x x x o o o o o ] /
 [         x x x x x x ]           1 2 3 4 5 6 7 8 9 0
   1 2 3 4 5 6 7 8 9 0

 \end{verbatim}
 */
class MatrixGenBanded	: public MatrixOp
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
  MatrixGenBanded(
    size_type                         m                       = 0
    ,size_type                        n                       = 0
    ,size_type                        kl                      = 0
    ,size_type                        ku                      = 0
    ,DMatrixSlice                   *MB                     = NULL
    ,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
    );

  // ///////////////////////////
    // Access representation

  /** \brief Initialize
   *
   * If called with all of the default arguments then #this# will become uninitialized.
   *
   * ToDo: Finish pre and post conditions!
   *
   * @param  m        [in] Determines the size of the banded matrix (m x n).  If
   *                  If #m == 0# then all of the following arguments should be left at
   * @param  n        [in] Determines the size of the banded matrix (m x n).
   * @param  kl       [in] Determines the lower band width of the matrix as defined by xGBMV(...).
   * @param  ku       [in] Determines the band width of the matrix as defined by xGBMV(...).
   * @param  MB       [in/state] If #MB != NULL# then this matrix (size (kl+ku+1) x n) is used to store
   *                  the original banded matrix #M# in the format of xGBMV(...).  This matrix must
   *                  be initialized on input.
   * @param  MB_release_resource_ptr
   *                  [in] Only significant if #MB != NULL#.  Points to a resource to
   *                  be released when #MB# is no longer needed.
   */
  void initialize(
    size_type                         m                       = 0
    ,size_type                        n                       = 0
    ,size_type                        kl                      = 0
    ,size_type                        ku                      = 0
    ,DMatrixSlice                   *MB                     = NULL
    ,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
    );

  /** \brief . */
  size_type kl() const;
  /** \brief . */
  size_type ku() const;
  /** \brief Get view of MB.
   */
  DMatrixSlice& MB();
  /** \brief . */
  const DMatrixSlice& MB() const;

  // /////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
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

private:
  
  // /////////////////////////////
  // Private data members

  size_type                       m_;
  size_type                       n_;
  size_type                       kl_;
  size_type                       ku_;
  DMatrixSlice                  MB_;
  release_resource_ptr_t          MB_release_resource_ptr_;

  // /////////////////////////////
  // Private member functions

  void assert_initialized() const;

}; // end class MatrixGenBanded

// ///////////////////////////////////////////////////////
// Inline members for MatrixGenBanded

inline
size_type MatrixGenBanded::kl() const
{
  return kl_;
}

inline
size_type MatrixGenBanded::ku() const
{
  return ku_;
}

inline
DMatrixSlice& MatrixGenBanded::MB()
{
  return MB_;
}

inline
const DMatrixSlice& MatrixGenBanded::MB() const
{
  return MB_;
}

} // end namespace ConstrainedOptPack

#endif // MATRIX_GEN_BANDED_H
