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

#ifndef ALAP_VECTOR_SPACE_Thyra_HPP
#define ALAP_VECTOR_SPACE_Thyra_HPP

#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Thyra_VectorSpaceBase.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>VectorSpace</tt> adapter subclass for <tt>Thyra::VectorSpaceBase<value_type> </tt>.
 *
 * Note that the default copy constructor and assignment operators are
 * allowed which yield in shallow copy, not deep copy.
 */
class VectorSpaceThyra : public VectorSpace {
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  VectorSpaceThyra();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceThyra(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
    ,const inner_prod_ptr_t                                                  &inner_prod    = Teuchos::null
    );
  /** \brief Initalize given a smart pointer to a <tt>Thyra::VetorSpace</tt> object.
   *
   * @param  thyra_vec_spc  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
   * @param  inner_prod    [in] Smart pointer to an inner product.  If <tt>inner_prod.get()==NULL</tt>
   *                       then a <tt>InnerProductThyra</tt> object will be used which will
   *                       point to this.
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec_spc.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == thyra_vec_spc.get()</tt>
   * <li>[<tt>inner_prod.get()!=NULL</tt>]
   *     <tt>this->inner_prod().get()==inner_prod.get()</tt>
   * <li>[<tt>inner_prod.get()==NULL</tt>]
   *     <tt>dynamic_cast<const InnerProductThyra*>(this->inner_prod().get()).thyra_vec_spc().get()==thyra_vec_spc.get()</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
    ,const inner_prod_ptr_t                                                  &inner_prod    = Teuchos::null
    );
  /** \brief Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > set_uninitialized();
  /** \brief Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * If <tt>this->thyra_vec_spc().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_vec_spc()</tt> object.
   */
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc() const;

  //@}

  /** @name Overridden from VectorSpace */
  //@{

  /** \brief . */
  space_ptr_t clone() const;
  /** \brief . */
  bool is_compatible(const VectorSpace& vec_spc ) const;
  /** \brief . */
  bool is_in_core() const;
  /** \brief . */
  index_type dim() const;
  /** \brief . */
  vec_mut_ptr_t create_member() const;
  /** \brief . */
  space_fcty_ptr_t small_vec_spc_fcty() const;
  /** \brief . */
  multi_vec_mut_ptr_t create_members(size_type num_vecs) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  const Thyra::VectorSpaceBase<value_type>                              *thyra_vector_space;
#else
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >  thyra_vec_spc_;
#endif

}; // end class VectorSpaceThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >&
VectorSpaceThyra::thyra_vec_spc() const
{
  return thyra_vec_spc_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_Thyra_HPP
