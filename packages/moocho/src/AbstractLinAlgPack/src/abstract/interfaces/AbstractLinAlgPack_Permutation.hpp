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

#ifndef ALAP_PERMUTATION_H
#define ALAP_PERMUTATION_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface to permutation matrices.
 *
 * A \c Permutation object is used to permute the elements within
 * a vector.  It is not a general linear operator since it does not
 * map between vector spaces.  It only permutes elements within the same
 * vector space.
 */
class Permutation {
public:

  /** \brief . */
  virtual ~Permutation() {}

  /** @name Vector space */
  //@{
  
  /** \brief Return a reference to a vector space object that this permutation is compatible with.
   */
  virtual const VectorSpace& space() const = 0;
  
  //@}

  /** @name Information */
  //@{

  /// Return the dimension of the permutation.
  virtual size_type dim() const = 0;

  /// Returns true if \c this is the identity permutation \a I.
  virtual bool is_identity() const = 0;

  /// Prints debug type of information
  virtual std::ostream& output(std::ostream& out) const = 0;

  //@}

  /** @name Vector permutations */
  //@{

  /** \brief Permute a vector <tt>op(P)*x -> y</tt>
   *
   * @param  P_trans  [in] <tt>op(P) = P</tt> for <tt>P_trans == BLAS_Cpp::no_trans</tt> or
   *                  <tt>op(P) = P'</tt> for <tt>P_trans == BLAS_Cpp::trans</tt>.
   * @param  x        [in] Vector.
   * @param  y        [out] Vector.
   *
   * Preconditions:<ul>
   * <li> <tt>y != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>x.space().is_compatible(this->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>y->space().is_compatible(this->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   */
  virtual void permute( 
    BLAS_Cpp::Transp    P_trans
    ,const Vector       &x
    ,VectorMutable      *y
    ) const = 0;

  /** \brief Permute a vector <tt>op(P)*y -> y</tt>
   *
   * @param  P_trans  [in] <tt>op(P) = P</tt> for <tt>P_trans == BLAS_Cpp::no_trans</tt> or
   *                  <tt>op(P) = P'</tt> for <tt>P_trans == BLAS_Cpp::trans</tt>.
   * @param  y        [in/out] Vector.
   *
   * Preconditions:<ul>
   * <li> <tt>y != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>y->space().is_compatible(this->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   */
  virtual void permute( 
    BLAS_Cpp::Transp    P_trans
    ,VectorMutable      *y
    ) const = 0;

  //@}

}; // end class Permutation

} // end namespace AbstractLinAlgPack

#endif // ALAP_PERMUTATION_H
