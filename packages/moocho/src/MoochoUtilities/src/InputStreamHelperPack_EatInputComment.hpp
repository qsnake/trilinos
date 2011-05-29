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

#ifndef EAT_INPUT_COMMENT_H
#define EAT_INPUT_COMMENT_H

#include "Moocho_ConfigDefs.hpp"

namespace InputStreamHelperPack {

/** @name namespace InputStreamHelperPack
 *
 * @memo Basic input stream helper functions.
 */
//@{

/** \brief Discards comment lines from an input stream.
 *
 * Call this function to discard text that includes comment lines and
 * lines with only newline chars.  Here comment lines are in the tradition
 * of Fortran in that the comment line must begin with the comment
 * identification character (user supplied) and ends with a newline char
 * ('\n').
 *
 * In particular this function will discard any white
 * space before a set of consecutive comment lines and
 * the comment lines and will leave the input steam
 * at the first non-comment line that does not start with '\n'.
 */
std::istream& eat_comment_lines(std::istream& is, char comment_identifier);

//@}

}	// end namespace InputStreamHelperPack

#endif	// EAT_INPUT_COMMENT_H
