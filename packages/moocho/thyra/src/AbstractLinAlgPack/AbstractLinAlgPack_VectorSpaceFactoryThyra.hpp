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

#ifndef ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP
#define ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP

#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>VectorSpaceFactory</tt> adapter subclass for <tt>Thyra::VectorSpaceBase</tt>.
 */
class VectorSpaceFactoryThyra : public VectorSpaceFactory {
public:

  /** @name Constructors / Initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  VectorSpaceFactoryThyra();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceFactoryThyra( const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> > &thyra_vec_spc_fcty );
  /** \brief Initalize given a smart pointer to a <tt>Thyra::VetorSpaceFactory</tt> object.
   *
   * @param  thyra_vec_spc_fcty  [in] Smart pointer to Thyra vector
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec_spc_fcty.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc_fcty().get() == thyra_vec_spc_fcty.get()</tt>
   * </ul>
   */
  void initialize( const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> > &thyra_vec_spc_fcty );
  /** \brief Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorSpaceBase</tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc_fcty().get() == NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> > set_uninitialized();
  /** \brief Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorSpaceBase</tt> object.
   *
   * If <tt>this->thyra_vec_spc_fcty().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_vec_spc_fcty()</tt> object.
   */
  const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> >& thyra_vec_spc_fcty() const;

  //@}

  /** @name Overridden from VectorSpaceFactory */
  //@{

  /** \brief . */
  space_ptr_t create_vec_spc(index_type dim) const;

  //@}
  
private:

#ifdef DOXYGEN_COMPILE
  const Thyra::VectorSpaceFactoryBase<value_type>                             *thyra_vector_space_factory;
#else
  Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> >  thyra_vec_spc_fcty_;
#endif

}; // end class VectorSpaceFactoryThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> >&
VectorSpaceFactoryThyra::thyra_vec_spc_fcty() const
{
  return thyra_vec_spc_fcty_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP
