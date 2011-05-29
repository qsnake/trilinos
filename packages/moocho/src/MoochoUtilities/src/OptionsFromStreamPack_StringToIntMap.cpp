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

#include "OptionsFromStreamPack_StringToIntMap.hpp"

namespace OptionsFromStreamPack {

StringToIntMap::StringToIntMap(  const std::string& name, int n, const char* strings[] )
  : name_(name)
{
  typedef map_t::value_type val_t;
  for( int i = 0; i < n; ++i ) {
    const bool unique = map_.insert( val_t( strings[i], i ) ).second;
    TEST_FOR_EXCEPTION(
      !unique,AlreadyExists
      ,"StringToIntMap::StringToIntMap(...): "
      << "Error, the option \"" << strings[i] << "\" is a duplicate for options_group \""
      << name_ << "\"" );
  }
}

int StringToIntMap::operator()( const std::string& str ) const
{
  map_t::const_iterator itr = map_.find( str );
  TEST_FOR_EXCEPTION(
    itr == map_.end(), DoesNotExist
    ,"StringToIntMap::operator(...): "
    << "Error, the option \"" << str << "\" is not recongnised for options_group \""
    << name_ << "\"" );
  return (*itr).second;	
}

} // end namespace OptionsFromStreamPack
