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
// These are subroutines that are called by QPKWIK to output data for
// debugging purposes.

#include <ostream>

namespace QPKWIK_Output {

/// 
/** Class for setting the output stream (destructor unsets it).
  *
  * Create an object of this type to set a stream for outputing
  * before you call QPKWIK.
  * This is not tread safe.
  */
class set_output {
public:
  /** \brief . */
  set_output(std::ostream* out);
  /** \brief . */
  ~set_output();
private:
  // not defined and not to be called
  set_output();
  set_output(const set_output&);
  set_output& operator=(const set_output&);
};	// end class set_output

// Output stream to use (default == 0, no output).
extern std::ostream* out;

}	// end namespace QPKWIK_Output
