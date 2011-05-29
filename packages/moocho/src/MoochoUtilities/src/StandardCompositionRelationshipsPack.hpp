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

#ifndef STANDARD_COMPOSITION_RELATIONSHIPS_PACK_H
#define STANDARD_COMPOSITION_RELATIONSHIPS_PACK_H

#include "Moocho_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"

namespace StandardCompositionRelationshipsPack {

/** @name <<std comp>> Stereotype Implementation Helper Package.
  *
  * This is the set of helper functions specified in the diagram
  * "Class Diagram : <<std comp>> Stereotype Implementation Helper Package".
  */
//@{

/// Thrown when the reference has not been set.
class NoRefSet : public std::logic_error
{public: NoRefSet(const std::string& what_arg) : std::logic_error(what_arg) {}};

// Throw a NoRefSet exception
inline void ThrowNoRefSet(const char func_name[], const char name[])
{
  TEST_FOR_EXCEPTION(
    true,NoRefSet
    ,func_name << ": The reference for \'" << name << "\' has not been set yet"
    );
}

/// Assert that the reference is set.
template<class ContainedClass>
inline void assert_role_name_set(const ContainedClass* role_name_, const char func_name[]
  , const char name[])
{
  if(!role_name_) ThrowNoRefSet(func_name, name);
}

/** \brief . */
template<class ContainedClass>
inline void set_role_name(ContainedClass*& role_name_, bool& owns_role_name_, const char name[]
  , ContainedClass* role_name, bool owns_role_name)
{
  if(owns_role_name_ && role_name_ != role_name) delete role_name_;
  role_name_ = role_name; owns_role_name_ = owns_role_name;
}

/** \brief . */
template<class ContainedClass>
inline ContainedClass* get_role_name(ContainedClass* role_name_, bool owns_role_name_
  , const char name[])
{
  return role_name_;
}

/** \brief . */
template<class ContainedClass>
inline void set_owns_role_name(ContainedClass*& role_name_, bool& owns_role_name_
  , const char name[], bool owns_role_name)
{
  assert_role_name_set(role_name_, "set_owns_role_name()", name);
  owns_role_name_ = owns_role_name;
}

/** \brief . */
template<class ContainedClass>
inline bool owns_role_name(ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "owns_role_name()", name);
  return owns_role_name_;
}

/** \brief . */
template<class ContainedClass>
inline ContainedClass& role_name(ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "role_name()", name);
  return *role_name_;
}

/** \brief . */
template<class ContainedClass>
inline const ContainedClass& role_name(const ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "role_name()", name);
  return *role_name_;
}

/** \brief . */
template<class ContainedClass>
inline const ContainedClass& const_role_name(const ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "role_name()", name);
  return *role_name_;
}


/** \brief . */
template<class ContainedClass>
inline void destory_container_obj(ContainedClass* role_name_, bool owns_role_name_)
{
  if(owns_role_name_) delete role_name_;
}

//@}

}	// end namespace StandardCompositionRelationshipsPack

#endif // STANDARD_COMPOSITION_RELATIONSHIPS_PACK_H
