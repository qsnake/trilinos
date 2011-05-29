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

#ifndef MATRIX_CONVERT_TO_SPARSE_ENCAP_H
#define MATRIX_CONVERT_TO_SPARSE_ENCAP_H

#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief Sparse conversion subclass based on views of a \c MatrixExtractSparseElements object.
 *
 * ToDo:Finish documentation!
 */
class MatrixConvertToSparseEncap
  : virtual public MatrixConvertToSparse
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<const MatrixExtractSparseElements>  mese_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const IVector>                      i_vector_ptr_t;

  //@}

  /** @name Constructors/initializers */
  //@{

  /** \brief Construct to uninitialized.
   */
  MatrixConvertToSparseEncap();

  /** \brief Calls \c this->initialize().
   */
  MatrixConvertToSparseEncap(
    const mese_ptr_t           &mese
    ,const i_vector_ptr_t      &inv_row_perm
    ,const Range1D             &row_rng
    ,const i_vector_ptr_t      &inv_col_perm
    ,const Range1D             &col_rng
    ,const BLAS_Cpp::Transp    mese_trans
    ,const value_type          alpha = 1.0
    );

  /** \brief Initialize a permuted view of a sparse matrix.
   *
   * <tt>A = alpha * op( (P'*B*Q)(row_rng,col_rng) )</tt>
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const mese_ptr_t           &mese
    ,const i_vector_ptr_t      &inv_row_perm
    ,const Range1D             &row_rng
    ,const i_vector_ptr_t      &inv_col_perm
    ,const Range1D             &col_rng
    ,const BLAS_Cpp::Transp    mese_trans
    ,const value_type          alpha = 1.0
    );

  /** \brief Set uninitialized.
   *
   * ToDo: Finish documentation!
   */
  void set_uninitialized();

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  const mese_ptr_t& mese() const;
  /** \brief . */
  const i_vector_ptr_t& inv_row_perm() const;
  /** \brief . */
  const Range1D& row_rng() const;
  /** \brief . */
  const i_vector_ptr_t& inv_col_perm() const;
  /** \brief . */
  const Range1D& col_rng() const;
  /** \brief . */
  const BLAS_Cpp::Transp mese_trans() const;
  /** \brief . */
  const value_type alpha() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;

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

private:

  typedef Teuchos::RCP<const VectorSpace> space_ptr_t;  

#ifdef DOXYGEN_COMPILE
  const MatrixExtractSparseElements    *mese;
  const DenseLinAlgPack::IVector            *inv_row_perm;
  Range1D                              row_rng;
  const DenseLinAlgPack::IVector            *inv_col_perm;
  Range1D                              col_rng;
#else
  mese_ptr_t          mese_;
  BLAS_Cpp::Transp    mese_trans_;
  value_type          alpha_;
  Range1D             row_rng_;
  Range1D             col_rng_;
  i_vector_ptr_t      inv_row_perm_;
  i_vector_ptr_t      inv_col_perm_;
  size_type           nz_full_;
  space_ptr_t         space_cols_;
  space_ptr_t         space_rows_;
#endif

};	// end class MatrixConvertToSparseEncap

// /////////////////////////////////////////////
// Inline members

// Access

inline
const MatrixConvertToSparseEncap::mese_ptr_t&
MatrixConvertToSparseEncap::mese() const
{
  return mese_;
}

inline
const MatrixConvertToSparseEncap::i_vector_ptr_t&
MatrixConvertToSparseEncap::inv_row_perm() const
{
  return inv_row_perm_;
}

inline
const Range1D& MatrixConvertToSparseEncap::row_rng() const
{
  return row_rng_;
}

inline
const MatrixConvertToSparseEncap::i_vector_ptr_t&
MatrixConvertToSparseEncap::inv_col_perm() const
{
  return inv_col_perm_;
}

inline
const Range1D& MatrixConvertToSparseEncap::col_rng() const
{
  return col_rng_;
}

inline
const BLAS_Cpp::Transp
MatrixConvertToSparseEncap::mese_trans() const
{
  return mese_trans_;
}

inline
const value_type MatrixConvertToSparseEncap::alpha() const
{
  return alpha_;
}
  
}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_CONVERT_TO_SPARSE_ENCAP_H
