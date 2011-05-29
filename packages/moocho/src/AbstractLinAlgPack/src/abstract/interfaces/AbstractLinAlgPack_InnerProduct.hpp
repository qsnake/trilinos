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

#ifndef ALAP_INNER_PRODUCT_H
#define ALAP_INNER_PRODUCT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface for inner products.
 *
 * ToDo: Finish documentaion
 */
class InnerProduct {
public:

  /** \brief . */
  virtual ~InnerProduct() {}

  /** \brief Compute the inner product of two vectors.
   *
   * Preconditions:<ul>
   * <li> ToDo: Spell out
   * </ul>
   *
   * Postconditions:<ul>
   * <li> ToDo: Spell out
   * </ul>
   *
   * @param  v1  [in] First vector
   * @param  v2  [in] Second vector
   *
   * @return  Returns some inner product of two vectors within a vector space..
   */
  virtual value_type inner_prod(const Vector& v1, const Vector& v2) const = 0;

}; // end class InnerProduct

} // end namespace AbstractLinAlgPack

#endif  // ALAP_INNER_PRODUCT_H
