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
#ifndef SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H
#define SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H

#include "AbstractLinAlgPack_MatrixSymOpSerial.hpp"
#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all serial symmetric diagonal matrices with
 * significant zeros along the diagonal.
 */
class MatrixSymDiagSparse
  : virtual public MatrixSymOpSerial
  , virtual public MatrixConvertToSparse
{
public:

  /** \brief <<std member comp>> members for how many updates to compute
    * at once in the operation M_MtMtM(....).
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_updates_at_once );

  /** \brief The default value of num_updates_at_once == 0 is set to allow
    * this class to determine the appropriate size internally.
    */
  MatrixSymDiagSparse();

  /** @name To be overridden by subclass */
  //@{

  /// Give access to the sparse diagonal
  virtual const SpVectorSlice diag() const = 0;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  std::ostream& output(std::ostream& out) const;

  //@}

  /** @name Overridden from MatrixOpSerial */
  //@{

  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;

  //@}

  /** @name Overridden from MatrixSymOpSerial */
  //@{

  /** \brief Computes the dense symmetric matrix B += a*op(A')*M*op(A).
   *
   * This matrix is computed using a set of rank-1 updates.
   *
   * Runtime ~ O( (m^2)*nz )
   *
   * Storage ~ O( num_updates_at_once * m )
   *
   * Where:<ul>
   * <li> <tt>n = A.rows() == this->rows()</tt>
   * <li> <tt>m = A.cols()</tt>
   * <li> <tt>nz = this->diag().nz()</tt>
   * </ul>
   *
   * Note that a necessary condition for \c B to be full rank is for
   * <tt>nz >= m</tt>.
   *
   * Also note that this default implementation is only for nonnegative
   * diagonal entries.
   */
  void Mp_StMtMtM( DMatrixSliceSym* sym_lhs, value_type alpha
    , EMatRhsPlaceHolder dummy_place_holder
    , const MatrixOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
    , value_type beta ) const;

  //@}

  /** @name Overridden from MatrixConvertToSparse */
  //@{

  /** \brief . */
  index_type num_nonzeros(
    EExtractRegion        extract_region
    ,EElementUniqueness   element_uniqueness
    ) const;
  /** \brief . */
  void coor_extract_nonzeros(
    EExtractRegion                extract_region
    ,EElementUniqueness           element_uniqueness
    ,const index_type             len_Aval
    ,value_type                   Aval[]
    ,const index_type             len_Aij
    ,index_type                   Arow[]
    ,index_type                   Acol[]
    ,const index_type             row_offset
    ,const index_type             col_offset
    ) const;

  //@}

};	// end class MatrixSymDiagSparse

}	// end namespace AbstractLinAlgPack

#endif	// SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H
