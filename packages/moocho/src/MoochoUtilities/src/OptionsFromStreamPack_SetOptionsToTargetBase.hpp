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

#ifndef SET_OPTIONS_TO_TARGET_BASE_H
#define SET_OPTIONS_TO_TARGET_BASE_H

#include "StandardCompositionRelationshipsPack.hpp"

namespace OptionsFromStreamPack {

/** \brief Templated node class manipulating a reference to a target object
  * who will have its options set..
  */
template< class T >
class SetOptionsToTargetBase {
public:

  /// 
  SetOptionsToTargetBase( T* target = 0 )
    : target_(target)
  {}

  /** @name <<std aggr>> stereotype members for target.
    */
  //@{

  /** \brief . */
  void set_target(T* target)
  {	target_ = target; }
  /** \brief . */
  T* get_target()
  {	return target_; }
  /** \brief . */
  const T* get_target() const
  {	return target_; }
  /** \brief . */
  T& target()
  {	return StandardCompositionRelationshipsPack::role_name(target_, false, "target"); }
  /** \brief . */
  const T& target() const
  {	return StandardCompositionRelationshipsPack::role_name(target_, false, "target"); }

  //@}

private:
  T* target_;

};	// end class SetOptionsToTargetBase

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_TO_TARGET_BASE_H
