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

#ifndef PRINT_SUB_VECTOR_H
#define PRINT_SUB_VECTOR_H

#include "RTOpPack_OldTypes.hpp"

namespace RTOpPack {

/** \brief Print an <tt>RTOp_SubVector</tt>.
 *
 * @param  o           [in/out] Stream that output is sent to.  On value of
 *                     <tt>o.width()</tt> on input is used to set the width
 *                     for the columns that the vector elements are printed
 *                     into.  However, no matter what the width argument is
 *                     at least one space is guarantteed between vector elements
 *                     Of course, the other stream properties will be
 *                     used also (i.e. precision, scientific etc.) when
 *                     printing the elements.
 * @param  v           [in] sub-vector object who's elements are printed.
 * @param  print_dim   [in] If \c true, then the dimension of the sub-vector will
 *                     be printed before the vector elements on the following line
 *                     (see below).
 * @param  newline     [in] If \c true, then a newline will be printed after the last
 *                     vector elements is printed.
 *
 * The indexes used for the output elements are exactly as is explained in the
 * specification for <tt>RTOp_SubVector</tt>.  Each of the vector elements are
 * printed in sparse format (even if the vector is dense) with the value followed
 * by a \c : followed by the index.  For example, suppose a sparse sub-vector is printed with
 * <tt>global_offset == 5</tt>, <tt>sub_dim == 8</tt> and <tt>sub_nz == 4</tt>.
 * With <tt>o.precision() == 4</tt> and <tt>o.width() == 16</tt>, <tt>print_dim == true</tt>
 * and <tt>newline == true</tt>, the output sent to \c o would look something like:
 \verbatim
 
 8
    -1.3457e+01:6     8.352e-01:8   -5.9245e-03:9   -3.658e+02:11 
 \endverbatim
 *
 * Above, if <tt>print_dim == false</tt>, then the sub-vector dimension <tt>"8"</tt> and the
 * <tt>"\\n"</tt> after it would not have been printed.  Also, if <tt>newline == false</tt> 
 * , then the newline <tt>"\\n"</tt> after <tt>-3.658e+02:11<\tt> would not have been printed.
 * ToDo: Finish documentation!
 */
std::ostream& output(
  std::ostream& o, const SubVector& v
  ,bool print_dim , bool newline
  );

} // end namespace RTOpPack

#endif // PRINT_SUB_VECTOR_H
