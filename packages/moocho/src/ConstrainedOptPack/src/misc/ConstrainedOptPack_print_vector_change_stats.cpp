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

#include "ConstrainedOptPack_print_vector_change_stats.hpp"
#include "ConstrainedOptPack_vector_change_stats.hpp"

void ConstrainedOptPack::print_vector_change_stats(
    const DVectorSlice& x, const char x_name[]
  , const DVectorSlice& d, const char d_name[], std::ostream& out )
{
  value_type	max_term,	min_term,	av_term;
  size_type	max_k,		min_k;
  vector_change_stats(
      x, d
    , &max_term, &max_k
    , &min_term, &min_k
    , &av_term	);
  out	<< "\nmax(|"<<d_name<<"(i)|/(1+|"<<x_name<<"(i)|)"
      << " => |"<<d_name<<"("<<max_k<<")|/(1+|"<<x_name<<"("<<max_k<<")| = "<< max_term
    << "\nmin(|"<<d_name<<"(i)|/(1+|"<<x_name<<"(i)|)"
      << " => |"<<d_name<<"("<<min_k<<")|/(1+|"<<x_name<<"("<<min_k<<")| = "<< min_term
    << "\naverage(|"<<d_name<<"(i)|/(1+|"<<x_name<<"(i)|) = " << av_term << std::endl;
}
