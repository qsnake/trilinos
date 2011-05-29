#if 0
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

#include "AbstractLinAlgPack_sparse_bounds.hpp"

/** \brief Count the number of sparse bounds where at least one bound is
  * finite.
  */
AbstractLinAlgPack::size_type
AbstractLinAlgPack::num_bounds( const SpVectorSlice& bl, const SpVectorSlice& bu )
{
  SpVectorSlice::const_iterator
    bl_itr			= bl.begin(),
    bl_itr_end		= bl.end(),
    bu_itr			= bu.begin(),
    bu_itr_end		= bu.end();
  size_type num_bounds = 0;
  while( bl_itr != bl_itr_end || bu_itr != bu_itr_end ) {
    if( ( bl_itr != bl_itr_end )
      && ( bu_itr == bu_itr_end || bl_itr->indice() + bl.offset() < bu_itr->indice() + bu.offset() ) )
    {
      // Only the lower bound is finite
      ++bl_itr;
    }
    else if( ( bu_itr != bu_itr_end )
      && ( bl_itr == bl_itr_end || bu_itr->indice() + bu.offset() < bl_itr->indice() + bl.offset()) )
    {
      // Only the upper bound is finite
      ++bu_itr;
    }
    else if(bl_itr->indice() == bu_itr->indice()) {
      // Both bounds exist.
      ++bl_itr;
      ++bu_itr; 
    }
    ++num_bounds;
  }
  return num_bounds;
}


#endif // 0
