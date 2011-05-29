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

#ifndef ALAP_MATRIX_ZERO_H
#define ALAP_MATRIX_ZERO_H

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Implementation of a matrix with all zeros.
 *
 * This may seem like a silly class but it is helpful in some circumstances.
 * This class needs to be updated whenever methods are added or removed
 * from \c MatrixOp.
 */
class MatrixZero : public MatrixOp {
public:

  /** @name Constructors/initializers */
  //@{

  /// Calls <tt>this->initalize()</tt>
  MatrixZero(
    const VectorSpace::space_ptr_t&     space_cols = Teuchos::null
    ,const VectorSpace::space_ptr_t&    space_rows = Teuchos::null
    );

  /** \brief Initialize (or initialize) given the columns and rows vector spaces.
   *
   * Preconditions:<ul>
   * <li> [<tt>space_cols.get() == NULL</tt>] <tt><tt>space_rows.get() == NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>space_cols.get() != NULL</tt>] <tt><tt>space_rows.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>space_cols.get() == NULL</tt>] <tt>this->rows() == 0</tt>
   * <li> [<tt>space_cols.get() == NULL</tt>] <tt>this->cols() == 0</tt>
   * <li> [<tt>space_cols.get() != NULL</tt>] <tt>&this->space_cols() == space_cols.get()</tt>
   * <li> [<tt>space_cols.get() != NULL</tt>] <tt>&this->space_rows() == space_rows.get()</tt>
   * </ul>
   */
  void initialize(
    const VectorSpace::space_ptr_t&    space_cols
    ,const VectorSpace::space_ptr_t&   space_rows
    );

  //@}

  /** @name Overridden from MatrixBase */
  //@{
  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;
  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief . */
  void zero_out();
  /** \brief . */
  void Mt_S(value_type alpha);
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& M);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  bool Mp_StM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs
    ) const;
  /** \brief . */
  bool Mp_StMtP(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp M_trans
    ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ) const;
  /** \brief . */
  bool Mp_StPtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ,BLAS_Cpp::Transp M_trans
    ) const;
  /** \brief . */
  bool Mp_StPtMtP(
    MatrixOp* mwo_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_trans
    ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
    ) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* vs_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_rhs2_trans
    ,const Vector& v_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* vs_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_rhs2_trans
    ,const SpVectorSlice& sv_rhs3, value_type beta) const;
  /** \brief . */
  value_type transVtMtV(
    const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
    ,const Vector& v_rhs3) const;
  /** \brief . */
  value_type transVtMtV(
    const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const;
  /** \brief . */
  void syr2k(
     BLAS_Cpp::Transp M_trans, value_type alpha
    ,const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
    ,const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
    ,value_type beta, MatrixSymOp* symwo_lhs ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;
  /** \brief . */
  bool syrk(
     BLAS_Cpp::Transp M_trans, value_type alpha
    ,value_type beta, MatrixSymOp* sym_lhs ) const;

  //@}

private:

  VectorSpace::space_ptr_t  space_cols_;
  VectorSpace::space_ptr_t  space_rows_;

  //
  void assert_initialized() const;

}; // end class MatrixZero

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_ZERO_H
