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

#ifndef RELEASE_RESOURCE_REF_COUNT_PTR_H
#define RELEASE_RESOURCE_REF_COUNT_PTR_H

#include "ReleaseResource.hpp"
#include "Teuchos_RCP.hpp"

namespace MemMngPack {

/** \brief Template class that implements ReleaseResource interface for
 * a RCP<T> object.
 *
 * Note that ~ReleaseResource_ref_count_ptr() does not need to be
 * implemented since the compiler generated version will already
 * be correct.
 */
template <class T>
class ReleaseResource_ref_count_ptr : public ReleaseResource {
public:

  /** \brief . */
  typedef Teuchos::RCP<T>   ptr_t;

  /// Just give public access to pointer
  ptr_t  ptr;

  /// Construct from a pointer
  ReleaseResource_ref_count_ptr(const ptr_t& ptr);

  // ////////////////////////////////////
  // Overriddend from ReleaseResource

  /** \brief . */
  bool resource_is_bound() const;

}; // end class ReleaseResource_ref_count_ptr

// ////////////////////////////////
// Inline function definitions

template <class T>
inline
ReleaseResource_ref_count_ptr<T>::ReleaseResource_ref_count_ptr(const ptr_t& p)
  : ptr(p)
{}

// ///////////////////////////////
// Template function definitions

template <class T>
bool ReleaseResource_ref_count_ptr<T>::resource_is_bound() const
{
  return ptr.get() != 0;
}

} // end namespace MemMngPack

#endif // RELEASE_RESOURCE_REF_COUNT_PTR_H
