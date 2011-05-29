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

#ifndef VECTOR_CHANGE_STATS_H
#define VECTOR_CHANGE_STATS_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Compute statistics for change in a vector.
  *
  * Given two vectors x and d where we wish to generate statistics
  * for the update x+d this function computes the following
  * quantitines:
  *
  * max( |d(i)|/(1+|x(i)|), i=1...n )
  *   => #max_k# = k, #max_term# = |d(k)|/(1+|x(k)|) <= 1\\
  * #min( |d(i)|/(1+|x(i)|), i=1...n ) => #min_k# = k, #min_term# = |d(k)|/(1+|x(k)|)#\\
  * #average( |d(i)|/(1+|x(i)|), i=1...10 )# => #av_term#\\
  * 
  * The purpose of generating these satistics is to determine
  * by how much x+d differs from x.
  *
  * If |d(i)|/|x(i)| < mach_eps with x(i) > 0 then we know that d(i) will
  * be lost when added to x(i) so x(i) + d(i) == x(i).
  *
  */
void vector_change_stats( const DVectorSlice& x, const DVectorSlice& d
  , value_type* max_term, size_type* max_k
  , value_type* min_term, size_type* min_k
  , value_type* av_term );

}	// end namespace ConstrainedOptPack

#endif	// VECTOR_CHANGE_STATS_H
