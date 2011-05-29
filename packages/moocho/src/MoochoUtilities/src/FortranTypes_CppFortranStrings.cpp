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
//
// ToDo: 8/9/99: Make sure that these are valid ASCII conversions
// on each target platform.  For example, UNIX does not use
// all 8 bits for characters.


#include "Moocho_Config.h"


#ifdef HAVE_MOOCHO_FORTRAN


#include "FortranTypes_CppFortranStrings.hpp"


int FortranTypes::convert_to_f_int_string( const char string[], f_int i_string[]
  , f_int* string_len )
{
  *string_len = 0;
  for( ; *string != 0; ++(*string_len) )
    *i_string++ = *string++;	// Sse standard C conversions by default
  return 0;	// success
}


int FortranTypes::convert_from_f_int_string( const f_int i_string[], f_int string_len
  , char string[] )
{
  for( f_int i = 0; i < string_len; )
    *string++ = *i_string++;	// Sse standard C conversions by default
  *string = 0; // Null terminate the target string.
  return 0;	// success
}


#endif // HAVE_MOOCHO_FORTRAN
