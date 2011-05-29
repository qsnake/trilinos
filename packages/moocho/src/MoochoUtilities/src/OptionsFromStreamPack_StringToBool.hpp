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

#ifndef STRING_TO_BOOL_H
#define STRING_TO_BOOL_H

#include "Moocho_ConfigDefs.hpp"

namespace OptionsFromStreamPack {

/** \brief Convert a string "true" or "false" into bool #true# or #false#.
  *
  * If the input string is not "true" or "false" then the exception
  * "InputException" will be thrown and the message will include
  * the name of the option this value is for.
  */
bool StringToBool( const char* opt_name, const char* str );

}	// end namespace OptionsFromStreamPack 

#endif	// STRING_TO_BOOL_H
