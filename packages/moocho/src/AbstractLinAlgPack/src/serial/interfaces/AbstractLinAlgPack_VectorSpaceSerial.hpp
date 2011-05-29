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

#ifndef VECTOR_SPACE_SERIAL_H
#define VECTOR_SPACE_SERIAL_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Subclass for serial vector space objects that create <tt>VectorMutableDense</tt>
 * vector and <tt>MultiVectorMutableDense</tt> multi-vector objects.
 *
 * The default constructor, copy constructor and assignment operators
 * are allowed since they have the correct semantics.
 */
class VectorSpaceSerial
  : public AbstractLinAlgPack::VectorSpace
{
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceSerial( size_type dim = 0 );

  /** \brief Initialize given the dimension of the vector space.
   *
   * @param  dim   [in] The dimension of the vector space.
   */
  void initialize( size_type dim );

  //@}

  /** @name Overridden from VectorSpece */
  //@{

  /** \brief Returns true if <tt>vec_space.dim() == this->dim()</tt>.
   *
   * The assumption here is that <tt>Vector::get_sub_vector()</tt>
   * and <tt>VectorMutable::get_sub_vector()</tt> can be used to implement
   * all of the methods on an SMP machine in an efficient manner.
   */
   bool is_compatible(const VectorSpace& vec_space) const;
  /// Returns true
  bool is_in_core() const;
  /// Returns 0 if uninitialized
  index_type dim() const;
  /// Returns a <tt>VectorSpaceFactorySerial</tt> object
  space_fcty_ptr_t small_vec_spc_fcty() const;
  /** \brief . */
  space_ptr_t clone() const;
  /// Returns a <tt>VectorMutableDense</tt> object.
  vec_mut_ptr_t create_member() const;
  /// Returns a <tt>MultiVectorMutableDense</tt> object.
  multi_vec_mut_ptr_t create_members(size_type num_vecs) const;
  /** \brief . */
  space_ptr_t sub_space(const Range1D& rng) const;
  /** \brief . */
  space_ptr_t space(
    const GenPermMatrixSlice  &P
    ,BLAS_Cpp::Transp         P_trans
    ) const;
  //@}

private:

  size_type     dim_;

}; // end class VectorSpaceSerial

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_SERIAL_H
