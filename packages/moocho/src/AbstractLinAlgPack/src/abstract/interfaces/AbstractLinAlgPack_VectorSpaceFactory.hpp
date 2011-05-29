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

#ifndef VECTOR_SPACE_FACTORY_H
#define VECTOR_SPACE_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface for objects that can create vector spaces of a specified dimension.
 *
 * ToDo: Finish documentation!
 */
class VectorSpaceFactory
{
public:

  /** \brief . */
  typedef Teuchos::RCP<const InnerProduct>   inner_prod_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const VectorSpace>    space_ptr_t;

  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  virtual ~VectorSpaceFactory();

  /// Calls \c inner_prod()
  VectorSpaceFactory( const inner_prod_ptr_t& inner_prod = Teuchos::null );

  /** \brief Initialize with an inner product object that will be given to vector.
   *
   * @param  inner_prod  [in] Smart pointer to inner product strategy object.
   *                     If <tt>inner_prod.get()==NULL</tt> then an
   *                     \c InnerProductDot object will be used instead.
   *
   * Postconditions:<ul>
   * <li> [<tt>inner_prod.get() != NULL</tt>] <tt>this->inner_prod().get() == inner_prod.get()</tt>
   * <li> [<tt>inner_prod.get() == NULL</tt>] <tt>dynamic_cast<InnerProductDot*>(this->inner_prod().get()) != NULL</tt>
   * </ul>
   */
  virtual void inner_prod( const inner_prod_ptr_t& inner_prod );

  /** \brief Return the smart pointer to the inner product strategy object.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * </ul>
   */
  virtual const inner_prod_ptr_t inner_prod() const;

  //@}

  /** @name Pure virtual functions that must be overridden */
  //@{

  /** \brief Create a vector space of the given dimension.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == dim</tt>
   * <li> [<tt>this->inner_prod().get() != NULL</tt>] <tt>this->inner_prod().get() == return->inner_prod().get()</tt>
   * <li> [<tt>this->inner_prod().get() == NULL</tt>] <tt>dynamic_cast<InnerProductDot*>(return->inner_prod().get()) != NULL</tt>
   * </ul>
   *
   * @return  Returns a smart reference counted pointer to a dynamically
   * allocated vector space object that can be used to create vector.
   */
  virtual space_ptr_t create_vec_spc(index_type dim) const = 0;

  //@}
  
private:
#ifdef DOXYGEN_COMPILE
  InnerProduct       *inner_prod;
#else
  inner_prod_ptr_t   inner_prod_;
#endif

}; // end class VectorSpaceFactory

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_SPACE_FACTORY_H
