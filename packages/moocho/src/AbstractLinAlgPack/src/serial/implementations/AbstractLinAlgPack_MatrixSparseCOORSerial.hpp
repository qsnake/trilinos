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

#ifndef MATRIX_SPARSE_COOR_SERIAL_H
#define MATRIX_SPARSE_COOR_SERIAL_H

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack_MatrixLoadSparseElements.hpp"
#include "AbstractLinAlgPack_MatrixExtractSparseElements.hpp"
#include "ReleaseResource.hpp"

namespace AbstractLinAlgPack {

/** \brief Coordinate matrix subclass.
 *
 * ToDo: Finish documentation!
 */
class MatrixSparseCOORSerial
  : virtual public AbstractLinAlgPack::MatrixOp
  , virtual public MatrixLoadSparseElements
  , virtual public MatrixExtractSparseElements
{
public:

  /** @name Public types */
  //@{
  
  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>     release_resource_ptr_t;
  /** \brief Subclass to delete dynamically allocated memory with \c delete[].
   *
   * This subclass can be used by the client to cause
   */
  class ReleaseValRowColArrays
    : public MemMngPack::ReleaseResource
  {
  public:
    /// Gives pointers to buffers to delete[]
    ReleaseValRowColArrays(
      value_type   *val
      ,index_type  *row_i
      ,index_type  *col_j
      )
      :val_(val)
      ,row_i_(row_i)
      ,col_j_(col_j)
      ,owns_mem_(true)
    {}
    // Calls delete on buffers if <tt>this->owns_memory() == true</tt>
    ~ReleaseValRowColArrays();
    /// Overridden from ReleaseResource
    bool resource_is_bound() const;
    /// Release ownership of memory
    void release_ownership() { owns_mem_ = false; }
    value_type* val()   { return val_; }
    index_type* row_i() { return row_i_; }
    index_type* col_j() { return col_j_; }
  private:
    value_type  *val_;
    index_type  *row_i_;
    index_type  *col_j_;
    bool        owns_mem_;
    // not defined and not to be called
    ReleaseValRowColArrays();
    ReleaseValRowColArrays(const ReleaseValRowColArrays&);
    ReleaseValRowColArrays& operator=(const ReleaseValRowColArrays&);
  }; // end class ReleaseValRowColArrays

  //@}

  /** @name Constructors / initializers */
  //@{

  /** \brief Let \c this allocate its own memory.
   */
  MatrixSparseCOORSerial();

  /** \brief Give memory to use to store nonzero elements.
   *
   * @param  max_nz  [in] Size of the buffers \c val[], \c row_i[] and \c col_j[].  It is allowed for
   *                 <tt>max_nz == 0</tt> in which case this matrix will not have any nonzero elements
   *                 (i.e. the zero matrix).
   * @param  val     [in] Buffer to store the nonzero entry values.  On input <tt>val[k], for k = 0...nz-1</tt>
   *                 must be set to valid nonzero entry values if <tt>rows > 0 && nz > 0</tt>.
   * @param  row_i   [in] Buffer to store the nonzero entry row indexes.  On input <tt>row_i[k], for k = 0...nz-1</tt>
   *                 must be set to valid nonzero entry row indexes if <tt>rows > 0 && nz > 0</tt>.
   * @param  col_j   [in] Buffer to store the nonzero entry row indexes.  On input <tt>col_j[k], for k = 0...nz-1</tt>
   *                 must be set to valid nonzero entry row indexes if <tt>rows > 0 && nz > 0</tt>.
   * @param  release_resource
   *                 [in] Smart pointer to a <tt>MemMngPack::ReleaseResouce</tt> object that will
   *                 deallocate the memory for the buffers \c val[], \c row_i[] and \c col_j[] when \c this
   *                 is finished with them.  It is allowed for <tt>release_resource.get() == NULL</tt> in which
   *                 case the client will be responsible to free this memory.  The client can use the
   *                 subclass \c ReleaseValRowColArrays to allow \c this to delete memory created with
   *                 \c new[] by setting <tt>release_resource = Teuchos::rcp(
   *                 new ReleaseValRowColArrays(val,row_i,col_j) )</tt> where \c val, \c row_i and \c col_j
   *                 where allocated with \c new[].  The \c release_resource object represents the sharing
   *                 of the data \c val[], \c row_i[] and \c col_j[].
   * @param  rows    [in] Number of rows in the matrix.  Default = 0 for fully uninitialized.
   * @param  cols    [in] Number of columns in the matrix.  Ignored if <tt>rows == 0</tt>.
   * @param  nz      [in] Number of nonzero elements in \c val[], \c row_i[] and \c col_j[]
   *                 already set.  Setting a value of <tt>nz > 0</tt> allows the client
   *                 to setup a matrix object with nonzero elements ready to go.  Ignored
   *                 if <tt>rows == 0</tt>.
   * @param  check_input
   *                 [in] If \c true and <tt>rows > 0 && nz > 0</tt> then the row and colunmn
   *                 indexes in \c row_i[] and \c col_j[] will be checked!  The defualt is \c false.
   *
   * Preconditions:<ul>
   * <li> <tt>max_nz > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>val != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>row_i != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>col_j != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>rows >= 0</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>rows > 0</tt>] <tt>cols > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>rows > 0</tt>] <tt>0 <= nz <= max_nz</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>nz > 0</tt>] <tt>1 <= row_i[k] <= rows, for k = 0...nz-1</tt> (throw std::invalid_argument)
   * <li> [<tt>nz > 0</tt>] <tt>1 <= col_j[k] <= cols, for k = 0...nz-1</tt> (throw std::invalid_argument)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->release_resource().get() == release_resource</tt>
   * <li> [<tt>rows > 0</tt>] <tt>this->rows() == rows</tt>
   * <li> [<tt>rows > 0</tt>] <tt>this->cols() == cols</tt>
   * <li> [<tt>rows > 0</tt>] <tt>this->nz() == nz</tt>
   * </ul>
   */
  void set_buffers(
    size_type                      max_nz
    ,value_type                    *val
    ,index_type                    *row_i
    ,index_type                    *col_j
    ,const release_resource_ptr_t  &release_resource
    ,size_type                     rows                = 0
    ,size_type                     cols                = 0
    ,size_type                     nz                  = 0
    ,bool                          check_input         = false
    );
    
  /** \brief Release all owned memory and make uninitialized.
   *
   * Postconditions:<ul>
   * <li> On the next call to <tt>this->reinitialize(rows,cols,max_nz)</tt>
   *      \c this will allocate its own memory.
   * <li> 
   * </ul>
   */
  void set_uninitialized();

  //@}

  /** @name Access */
  //@{

  /** \brief . */
  value_type* val();
  /** \brief . */
  const value_type* val() const;
  /** \brief . */
  index_type* row_i();
  /** \brief . */
  const index_type* row_i() const;
  /** \brief . */
  index_type* col_j();
  /** \brief . */
  const index_type* col_j() const;
  /** \brief . */
  const release_resource_ptr_t& release_resource() const;

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
  MatrixOp& operator=(const MatrixOp& M);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2, value_type beta) const;

  // ToDo: Add more method overrides as they are needed!
  
  //@}
  
  /** @name Overridden from MatrixLoadSparseElements */
  //@{

  /** \brief . */
  void reinitialize(
    size_type                  rows
    ,size_type                 cols
    ,size_type                 max_nz
    ,EAssumeElementUniqueness  element_uniqueness
    );
  void reset_to_load_values();
  /** \brief . */
  void get_load_nonzeros_buffers(
    size_type      max_nz_load
    ,value_type    **val
    ,index_type    **row_i
    ,index_type    **col_j
    );
  /** \brief . */
  void commit_load_nonzeros_buffers(
    size_type      nz_commit
    ,value_type    **val
    ,index_type    **row_i
    ,index_type    **col_j
    );
  /** \brief . */
  void finish_construction( bool test_setup );

  //@}

  /** @name Overridden from MatrixExtractSparseElements */
  //@{

  /** \brief . */
  index_type count_nonzeros(
    EElementUniqueness    element_uniqueness
    ,const index_type     inv_row_perm[]
    ,const index_type     inv_col_perm[]
    ,const Range1D        &row_rng
    ,const Range1D        &col_rng
    ,index_type           dl
    ,index_type           du
    ) const;
  /** \brief . */
  void coor_extract_nonzeros(
    EElementUniqueness    element_uniqueness
    ,const index_type     inv_row_perm[]
    ,const index_type     inv_col_perm[]
    ,const Range1D        &row_rng
    ,const Range1D        &col_rng
    ,index_type           dl
    ,index_type           du
    ,value_type           alpha
    ,const index_type     len_Aval
    ,value_type           Aval[]
    ,const index_type     len_Aij
    ,index_type           Arow[]
    ,index_type           Acol[]
    ,const index_type     row_offset
    ,const index_type     col_offset
    ) const;

  //@}

private:

  // //////////////////////////////
  // Private types

  // //////////////////////////////
  // Public types

  size_type                 rows_;
  size_type                 cols_;
  size_type                 max_nz_;
  EAssumeElementUniqueness  element_uniqueness_;
  size_type                 nz_;
  value_type                *val_;
  index_type                *row_i_;
  index_type                *col_j_;
  release_resource_ptr_t    release_resource_;

  bool                      self_allocate_; // True if this allocates the memory

  VectorSpaceSerial         space_cols_;
  VectorSpaceSerial         space_rows_;

  size_type                 max_nz_load_;     // cashed
  bool                      reload_val_only_; // cashed
  size_type                 reload_val_only_nz_last_; // cashed

  // //////////////////////////////
  // Private member functions

  void make_storage_unique();

  // static
  static release_resource_ptr_t  release_resource_null_;

}; // end class MatrixSparseCOORSerial

// //////////////////////////////////
// Inline members

inline
value_type* MatrixSparseCOORSerial::val()
{
  make_storage_unique();
  return val_;
}

inline
const value_type* MatrixSparseCOORSerial::val() const
{
  return val_;
}

inline
index_type* MatrixSparseCOORSerial::row_i()
{
  make_storage_unique();
  return row_i_;
}

inline
const index_type* MatrixSparseCOORSerial::row_i() const
{
  return row_i_;
}

inline
index_type* MatrixSparseCOORSerial::col_j()
{
  make_storage_unique();
  return col_j_;
}

inline
const index_type* MatrixSparseCOORSerial::col_j() const
{
  return col_j_;
}

inline
const MatrixSparseCOORSerial::release_resource_ptr_t&
MatrixSparseCOORSerial::release_resource() const
{
  return self_allocate_ ? release_resource_null_ : release_resource_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SPARSE_COOR_SERIAL_H

