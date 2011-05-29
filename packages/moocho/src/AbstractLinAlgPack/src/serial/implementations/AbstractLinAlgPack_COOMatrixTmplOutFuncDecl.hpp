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

#ifndef COO_MATRIX_TMPL_OUT_FUNC_DECL_H
#define COO_MATRIX_TMPL_OUT_FUNC_DECL_H

#include "SparseLinAlgPackIOBasic.hpp"

namespace AbstractLinAlgPack {

/** \brief @name COO matrix output stream function.
  *
  * This is a functions that are used to output a templated COO matrix object
  * to a char based output stream.  The COOMatrixTemplateInterface specification
  * is used.
  *
  * The output format is diferent depending on whether the
  * bits #SparseLinAlgPackIO::ignore_dim_bit#,  #LinAlgPackIO::no_insert_newlines_bit#
  * , and #SparseLinAlgPackIO::ignore_nz_bit# are set.
  * The default output format is:
  *
  *	#rows  cols  nz#\\
  * #	val1:i1:j1	val2:i2:j2 ... valnz:inz:jnz#\\
  *
  * Each of the elements (val:i:j) are are put into columns according to the width set in the 
  * output stream #os# and other formating commands when it is called.  Even if the set
  * width is 0 or less than the number of char's for the element a space ' ' will be inserted
  * between them.  The elements are formated according to the format in the stream #os#.
  *
  * If #exta_flags & SparseLinAlgPackIO::ignore_dim_bit# == 0# then #rows# and #cols# will not
  * be output.
  *
  * If #exta_flags & SparseLinAlgPackIO::ignore_nz_bit# == 0# then #nz# will not
  * be output.
  *
  * If #exta_flags & LinAlgPackIO::no_insert_newlines_bit == 0# then a newline charactor
  * will not be inserted after the last element, otherwise (by default) it will be.
  *
  * If #coom.nz() == 0# then no elements will be output.
  *
  * If any of the output operations fails then a #std::ios_base::failure# exception is thrown. 
  */
template <class T_COOM>
std::ostream& output_COOM(std::ostream& os, const T_COOM& coom
  , SparseLinAlgPackIO::fmtflags extra_flags);

}	// end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_OUT_FUNC_DECL_H
