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

#ifndef GENMATRIX_IN_FUNC_H
#define GENMATRIX_IN_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DMatrix/DMatrixSlice input stream functions.
  *
  * These are functions that are used to read a DMatrix or DMatrixSlice object in from a
  * formated input stream.
  *
  * The input format is diferent depending on the on whether the bit
  * #LinAlgPackIO::ignore_dim_bit# is set.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0# then the input format is:
  *
  * Case 1\\
  *	#m  n#\\
  * #gm(1,1) gm(1,2) gm(1,3) ... gm(1,n)#\\
  * #gm(2,1) gm(2,2) gm(2,3) ... gm(2,n)#\\
  * #   .       .       .           .#\\
  * #gm(m,1) gm(m,2) gm(m,3) ... gm(m,n)#\\
  *
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the input format is:
  *
  * Case 2\\
  * #gm(1,1) gm(1,2) gm(1,3) ... gm(1,n)#\\
  * #gm(2,1) gm(2,2) gm(2,3) ... gm(2,n)#\\
  * #   .       .       .           .#\\
  * #gm(m,1) gm(m,2) gm(m,3) ... gm(m,n)#\\
  *
  * The numbers of the input must be seperated by white space (the line breaks are
  * for looks only) and be valid C language numeric constants.
  *
  * In addition, comment lines may be inserted between rows of the input matrix.
  * These comment lines take the form of Fortan comment lines in that they start on
  * new lines (after a '\n' char) with a '*' char and end at the end of a line
  * ('\n' terminated).  After the elements for a row are read in the function
  * #eat_comment_lines(is,'*');# is called.
  *
  * For example, the input format for the matrix {1.1 1.2; 2.1 2.2}
  * with comments for case 1 is:
  *
  * #2  2#\\
  * #* This is the first row#\\
  * #1.1	1.2#\\
  * #* This is the second row#\\
  * #2.1	2.1#\\
  *
  * And for case 2 is:
  *
  * #* This is the first row#\\
  * #1.1	1.2#\\
  * #* This is the second row#\\
  * #2.1	2.1#\\
  *
  * It is permisible for the dimension #m# and #n# in case 1 to be 0.
  * In this case there will be no elements.  So to input an empty matrix you would use:
  *
  * #0  0#\\
  *
  * If one of the dimenstions is zero but not the other then a #std::length_error# 
  * exception will be thrown.
  * If any of the input operations fails then a LinAlgPackIO::InputException exception
  * is thrown.  In other words if #is.fail()# or #is.eof()# is true
  * before all of the elements have been read in then the exception is thrown.
  * Also if the stream becomes corrupted (#is.bad() == true#) then a #std::ios_base::failure#
  * exception is thrown. 
  */
// @{

/** \brief . */
/* * DMatrix input stream function.
  *
  * Inputs a DMatrix object from an input stream.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0# then #gm# is resized to #m# x #n#
  * given in the file.  If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0#
  * then the number of elements read in depends on the current dimension of #gm#.
  */
std::istream& input(std::istream& is, DMatrix* gm, LinAlgPackIO::fmtflags extra_flags);

/** \brief . */
/* * DMatrixSlice input stream function.
  *
  * Inputs a DMatrixSlice object from an input stream.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0# then the dimension (sized) of #gms#
  * is compared to the #m# and #n# given in the file and if they are not equal
  * then a #LinAlgPackIO::InputException# is thrown.  If #gms# is unsized then it is resized
  * to #m# x #n#.  If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the number of
  * elements read in depends on the current size of #gms#.
  */
std::istream& input(std::istream& is, DMatrixSlice* gms, LinAlgPackIO::fmtflags extra_flags);

// @}

}	// end namespace DenseLinAlgPack

#endif	// GENMATRIX_IN_FUNC_H
