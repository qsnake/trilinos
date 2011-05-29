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

#ifndef VECTOR_OUT_FUNC_H
#define VECTOR_OUT_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DVectorSlice output stream function.
  *
  * This is a functions that are used to output a DVectorSlice object
  * to a char based output stream.
  *
  * The output format is diferent depending on the on whether the
  * bits #LinAlgPackIO::ignore_dim_bit# and #LinAlgPackIO::no_insert_newlines_bit# are set.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then the output format is:
  *
  * Case 1\\
  *	#v.size()#\\
  * #	v(1)	v(2)	v(3) ... v(v.size())#\\
  *
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the output format is:
  *
  * Case 2\\
  * #	v(1)	v(2)	v(3) ... v(v.size())#\\
  *
  * Each of the elements are are put into columns according to the width set in the 
  * output stream #os# and other formating commands when it is called.  Even if the set
  * width is 0 or less than
  * the number of char's for the element a space ' ' will be inserted between them.
  * The elements are formated according to the format in the stream #os#.
  *
  * If #exta_flags & LinAlgPackIO::no_insert_newlines_bit == 0# then a newline charactor
  * will not be inserted after the last element, otherwise (by default) it will be.
  *
  * If #vs.size() == 0# then no elements will be output and if
  * #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then only the operation
  * #os << vs.size() << endl;# will be performed.
  *
  * If any of the output operations fails then a #std::ios_base::failure# exception is thrown. 
  */
std::ostream& output(std::ostream& os, const DVectorSlice& vs, LinAlgPackIO::fmtflags extra_flags);

}	// end namespace DenseLinAlgPack

#endif	// VECTOR_OUT_FUNC_H
