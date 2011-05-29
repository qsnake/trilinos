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

#include <limits>

#include "ConstrainedOptPack_vector_change_stats.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_ScalarTraits.hpp"

void ConstrainedOptPack::vector_change_stats(
    const DVectorSlice& x, const DVectorSlice& d
  , value_type* max_term, size_type* max_k
  , value_type* min_term, size_type* min_k
  , value_type* av_term )
{
  typedef Teuchos::ScalarTraits<value_type> ST;
  DenseLinAlgPack::VopV_assert_sizes( x.dim(), d.dim() );
  const value_type
    min_num		= std::numeric_limits<value_type>::min(),
    inf			= std::numeric_limits<value_type>::max();
  // Initialize statistics
  *max_term	= 0.0;
  *max_k		= 0;
  *min_term	= inf;
  *min_k		= 0;
  *av_term	= 0.0;
  // Compute statistics
  DVectorSlice::const_iterator
    x_itr	= x.begin(),
    d_itr	= d.begin();
  for( size_type i = 1; x_itr != x.end(); ++i, ++d_itr, ++x_itr ) {
    // Compute this ratio and make sure we are not dividing by zero.
    // We only care about ratios less than 1.0 and also deal
    // with the case that both x(i) and d(i) are zero (in which
    // case the ratio should be zero since x(i) == x(i) + d(i)).
    const value_type
      term = ST::magnitude(*d_itr) / ( 1 + ST::magnitude(*x_itr) );
    if( term > *max_term ) {
      *max_term	= term;
      *max_k		= i;
    }
    if( term < *min_term ) {
      *min_term	= term;
      *min_k		= i;
    }
    *av_term += term;
  }
  *av_term = *av_term / x.dim();
}
