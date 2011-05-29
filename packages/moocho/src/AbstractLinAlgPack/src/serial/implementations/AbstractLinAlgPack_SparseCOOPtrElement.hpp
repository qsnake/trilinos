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

#ifndef SPARSE_COO_PTR_ELEMENT_H
#define SPARSE_COO_PTR_ELEMENT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Sparse pointer element type for a COO matrix (val, ivect, jvect).
 *
 * This class abstracts a sparse element of a templated
 * type from a coordinate matrix. It
 * has a pointer to the value of the element.
 *
 * The default assignment operator and copy constructor
 * are allowed.
 */
template <class T_Index, class T_Value>
class SparseCOOPtrElement {
public:
  /** @name Public Typedefs. */
  //@{

  /** \brief . */
  typedef T_Value						value_type;
  /** \brief . */
  typedef T_Index						index_type;

  //@}

  /** @name Constructors */
  //@{

  /// Construct uninitialized (poiner to value set to zero) (#index() == 0#).
  SparseCOOPtrElement() : pvalue_(0), row_i_(0), col_j_(0)
  {}

  /// Construct with a pointer to the value and index set
  SparseCOOPtrElement(value_type* pvalue, index_type row_i, index_type col_j)
    : pvalue_(pvalue), row_i_(row_i), col_j_(col_j)
  {}

  /// Initialize
  void initialize(value_type* pvalue, index_type row_i, index_type col_j) {
    pvalue_	= pvalue;
    row_i_	= row_i;
    col_j_	= col_j;
  }
  
  //@}

  /** @name Value and index access */
  //@{ 

  /** \brief . */
  value_type& value()
  {
    return *pvalue_;
  }
  /** \brief . */
  value_type value() const
  {
    return *pvalue_;
  }
  /** \brief . */
  index_type row_i() const
  {
    return row_i_;
  }
  /** \brief . */
  index_type col_j() const
  {
    return col_j_;
  }
  /// Change the indexs
  void change_indexes(index_type row_i, index_type col_j)
  {
    row_i_ = row_i;
    col_j_ = col_j;
  }
  /// Change the element pointer
  void change_value_ptr(value_type* pvalue)
  {
    pvalue_ = pvalue;
  }

  //@}
private:
  value_type*				pvalue_;
  index_type				row_i_, col_j_;

};	// end class SparseCOOPtrElement

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_COO_PTR_ELEMENT_H
