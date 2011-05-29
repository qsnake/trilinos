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

#ifndef SPARSE_PTR_ELEMENT_H
#define SPARSE_PTR_ELEMENT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Sparse pointer element type.
  *
  * This class abstracts a sparse element of a templated
  * type.  It is ment to be used in a sparse vector.  It
  * has a pointer to the value of the element.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  */
template <class T_Indice, class T_Value>
class SparsePtrElement {
public:
  /** @name Public Typedefs. */
  //@{

  /** \brief . */
  typedef T_Value							value_type;
  /** \brief . */
  typedef T_Indice						indice_type;

  //@}

  /** @name Constructors */
  //@{

  /// Construct uninitialized (poiner to value set to zero) (#indice() == 0#).
  SparsePtrElement() : indice_(0), pvalue_(0)
  {}

  /// Construct with a pointer to the value and indice set
  SparsePtrElement(indice_type indice, value_type* pvalue) : indice_(indice), pvalue_(pvalue)
  {}
  
  //@}

  /** @name Value and indice access */
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
  indice_type indice() const
  {
    return indice_;
  }
  /// Change the indice
  void change_indice(indice_type indice)
  {
    indice_ = indice;
  }
  /// Change the element pointer
  void change_value_ptr(value_type* pvalue)
  {
    pvalue_ = pvalue;
  }

  //@}
private:
  indice_type				indice_;
  value_type*				pvalue_;

};	// end class SparsePtrElement

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_PTR_ELEMENT_H
