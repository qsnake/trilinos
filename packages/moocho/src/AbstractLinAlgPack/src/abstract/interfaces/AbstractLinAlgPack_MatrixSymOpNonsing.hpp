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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_SYM_WITH_OP_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_SYM_WITH_OP_NONSINGULAR_H

#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixSymNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all polymorphic symmetrix nonsingular matrices that
 * can be used to compute matrix-vector products and solve for
 * linear systems relatively efficently.
 */
class MatrixSymOpNonsing 
  : virtual public MatrixSymOp
  , virtual public MatrixSymNonsing
  , virtual public MatrixOpNonsing
{
public:

  /** @name Public types */
  //@{

#ifndef DOXYGEN_COMPILE
  /** \brief . */
  typedef Teuchos::RCP<const MatrixSymOpNonsing>    mat_mswons_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<MatrixSymOpNonsing>          mat_mswons_mut_ptr_t;
#endif

  //@}

  /** @name Clone */
  //@{

  /** \brief Clone the non-const matrix object (if supported).
   *
   * The default implementation returns NULL which is perfectly acceptable.
   * A matrix object is not required to return a non-NULL value but almost
   * every good matrix implementation will.
   */
  virtual mat_mswons_mut_ptr_t clone_mswons();

  /** \brief Clone the const matrix object (if supported).
   *
   * The behavior of this method is the same as for the non-const version
   * above except it returns a smart pointer to a const matrix object.
   */
  virtual mat_mswons_ptr_t clone_mswons() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mut_ptr_t clone();
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_ptr_t clone() const;
  //@}

  /** @name Overridden from MatrixNonsing */
  //@{
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mns_mut_ptr_t clone_mns();
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mns_ptr_t clone_mns() const;
  //@}

  /** @name Overridden from MatrixSymOp */
  //@{
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mswo_mut_ptr_t clone_mswo();
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mswo_ptr_t clone_mswo() const;
  //@}

  /** @name Overridden from MatrixSymNonsing */
  //@{
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_msns_mut_ptr_t clone_msns();
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_msns_ptr_t clone_msns() const;
  //@}

  /** @name Overridden from MatrixOpNonsing */
  //@{
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mwons_mut_ptr_t clone_mwons();
  /// Returns <tt>this->clone_mswons()</tt>.
  mat_mwons_ptr_t clone_mwons() const;
  //@}

  /// Calls operator=(MatrixOp&)
  MatrixSymOpNonsing& operator=(const MatrixSymOpNonsing& M)
  { static_cast<MatrixOp*>(this)->operator=(M); return *this; }

}; // end class MatrixSymOpNonsing

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_SYM_WITH_OP_NONSINGULAR_H
