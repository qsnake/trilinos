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

#ifndef ALAP_VECTOR_SUB_VIEW_H
#define ALAP_VECTOR_SUB_VIEW_H

#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorApplyOpSerialBase.hpp"
#include "AbstractLinAlgPack_VectorSpaceSubSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Concrete subclass for a default sub-view implementation for a Vector
 * object.
 *
 * Not all of the methods from Vector are overridden, only those that
 * need to be or may result in better performance.
 *
 * There is really not much to this vector subclass.  The subclass is
 * only possible because of the <tt>first_ele</tt>, <tt>sub_dim</tt>,
 * and <tt>global_offset</tt> options with <tt>apply_op()</tt>.  The
 * vector space object returned by <tt>this->space()</tt> is of type
 * <tt>VectorSpaceSubSpace</tt> which in turn relys on <tt>VectorSpace::sub_space()</tt>.
 *
 * The default constructor and copy constructors are allowed but the default
 * assignment operator is not allowed.
 */
class VectorSubView
  : virtual public Vector
  , virtual protected VectorApplyOpSerialBase
{
public:

  /** \brief Constructs to uninitialized.
   *
   * Postconditions: see \c set_uninitialized().
   */
  VectorSubView();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSubView( const vec_ptr_t& full_vec, const Range1D& rng );
  /** \brief Initialize a sub-view based on a full vector.
   *
   * Constructs a view of the vector <tt>this = vec(rng)</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>full_vec.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>rng.full_range() == false</tt>] <tt>rng.lbound() <= full_vec->dim()</tt> (throw <tt>std::out_of_range</tt>).
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_ele(i) == full_vec->get_ele(rng.lbound()-1+i)</tt>,
   *      for <tt>i = 1...rng.size()</tt>
   * </ul>
   *
   * @param  full_vec  [in] The original full vector.
   * @param  rng       [in] The range of elements in <tt>full_vec</tt> that <tt>this</tt> vector will represent.
   */
  void initialize( const vec_ptr_t& full_vec, const Range1D& rng );
  /** \brief Set uninitialized()
   *
   * Postconditions:<ul>
   * <li> <tt>this->dim() == 0</tt>
   * <li> <tt>this->full_vec() = NULL</tt>
   * </ul>
   */
  void set_uninitialized();
  /** \brief . */
  const vec_ptr_t& full_vec() const;
  /** \brief . */
  const VectorSpaceSubSpace& space_impl() const;

  /** @name Overridden from Vector */
  //@{

  /** \brief . */
  const VectorSpace& space() const;
  /** \brief . */
  index_type dim() const;
  /** \brief Calls \c apply_op() on the underlying full vectors.
   *
   * Preconditions:<ul>
   * <li> <tt>dynamic_cast<const VectorSubView*>(vecs[k]) != NULL</tt>, for <tt>k=0..num_vecs</tt>
   *      (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>dynamic_cast<VectorMutableSubView*>(targ_vecs[k]) != NULL</tt>, for <tt>k=0..num_targ_vecs</tt>
   *      (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>dynamic_cast<const VectorSubView*>(vecs[k])->full_vec()->space().is_compatible(
   *      this->full_vec()->space() ) == true</tt>, for <tt>k=0..num_vecs</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>dynamic_cast<VectorMutableSubView>(targ_vecs[k])->full_vec()->space().is_compatible(
   *      this->full_vec()->space() ) == true</tt>, for <tt>k=0..num_targ_vecs</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   */
  void apply_op(
    const RTOpPack::RTOp& op
    ,const size_t num_vecs, const Vector* vecs[]
    ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
    ,RTOpPack::ReductTarget *reduct_obj
    ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
    ) const;
  /** \brief . */
  value_type get_ele(index_type i) const;
  /** \brief . */
  vec_ptr_t sub_view( const Range1D& rng ) const;
  /** \brief . */
  void get_sub_vector( const Range1D& rng, RTOpPack::SubVector* sub_vec ) const;
  /** \brief . */
  void free_sub_vector( RTOpPack::SubVector* sub_vec ) const;

  //@}

private:

  vec_ptr_t                  full_vec_;   ///< If full_vec_.get() == NULL, the vector is uninitalized (dim == 0).
  VectorSpaceSubSpace        space_;      ///< The space that this vector belongs to.

  // Not defined and not to be called
  VectorSubView& operator=(const VectorSubView&);
  
}; // end class VectorSubView

// /////////////////////////////////////////////
// Inline members

inline
VectorSubView::VectorSubView()
{}

inline
const VectorSubView::vec_ptr_t&
VectorSubView::full_vec() const
{
  return full_vec_;
}

inline
const VectorSpaceSubSpace& VectorSubView::space_impl() const
{
  return space_;
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_VECTOR_SUB_VIEW_H
