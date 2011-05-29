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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_RCP.hpp"

namespace OptionsFromStreamPack {
  class OptionsFromStream;
}

namespace AbstractLinAlgPack {

/** \brief Interface for a factory object that will create <tt>BasisSystem</tt> objects.
 *
 * 
 */
class BasisSystemFactory : public Teuchos::AbstractFactory<BasisSystem>
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const OptionsFromStreamPack::OptionsFromStream>             options_ptr_t;

  //@}

  /** \brief . */
  virtual ~BasisSystemFactory() {}

  /** \brief Set the options that will be used to determine what basis system will be returned
   * from <tt>this->create()</tt>.
   *
   * Note that it is allowed for the client to alter <tt>*options.get()</tt> after
   * this method is called so <tt>this</tt> had better read the options inside of the
   * <tt>this->create()</tt> method.
   */
  virtual void set_options( const options_ptr_t& options ) = 0;

  /** \brief Get the <tt>OptionsFromStream</tt> object being used to extract the options from.
   */
  virtual const options_ptr_t& get_options() const = 0;

}; // end class BasisSystemFactory

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H
