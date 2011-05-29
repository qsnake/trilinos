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

#ifndef ALAP_INNER_PRODUCT_Thyra_H
#define ALAP_INNER_PRODUCT_Thyra_H

#include "AbstractLinAlgPack_InnerProduct.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace AbstractLinAlgPack {


/** \brief Implements the inner product using
 * <tt>Thyra::VectorSpaceBase::scalarProd()</tt>.
 */
class InnerProductThyra : public InnerProduct {
public:

  /** @name Constructors / Initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  InnerProductThyra();

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  InnerProductThyra( const RCP<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc );

  /** \brief Initalize given a smart pointer to a <tt>Thyra::VetorSpace</tt> object.
   *
   * \param thyra_vec_spc [in] Smart pointer to Thyra vector
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec_spc.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == thyra_vec_spc.get()</tt>
   * </ul>
   */
  void initialize( const RCP<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc );

  /** \brief Set to uninitialized and return smart pointer to the internal
   * <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == NULL</tt>
   * </ul>
   */
  RCP<const Thyra::VectorSpaceBase<value_type> > set_uninitialized();

  /** \brief Return a (converted) smart pointer to the internal smart pointer
   * to the <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * If <tt>this->thyra_vec_spc().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_vec_spc()</tt> object.
   */
  const RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc() const;

  //@}

  /** @name Overridden from InnerProduct */
  //@{

  /** \brief . */
  value_type inner_prod(const Vector& v1, const Vector& v2) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  const Thyra::VectorSpaceBase<value_type> *thyra_vector_space;
#else
  RCP<const Thyra::VectorSpaceBase<value_type> > thyra_vec_spc_;
#endif

}; // end class InnerProductThyra


// ///////////////////////////////
// Inline functions


inline
const RCP<const Thyra::VectorSpaceBase<value_type> >&
InnerProductThyra::thyra_vec_spc() const
{
  return thyra_vec_spc_;
}


} // end namespace AbstractLinAlgPack


#endif  // ALAP_INNER_PRODUCT_Thyra_H
