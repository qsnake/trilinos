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

#include "AbstractLinAlgPack_sparse_bounds_diff.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

void AbstractLinAlgPack::imp_sparse_bnd_diff(
    int						sign
  , const SpVectorSlice		&sv
  , BLAS_Cpp::Uplo			uplo
  , const DVectorSlice			&v
  , DVectorSlice				*r
  )
{
  DenseLinAlgPack::Vp_V_assert_sizes(r->size(),sv.size());
  DenseLinAlgPack::VopV_assert_sizes(sv.size(),v.size());

  typedef DenseLinAlgPack::value_type value_type;
  const value_type
    inf = std::numeric_limits<value_type>::max();
  *r = ( uplo == BLAS_Cpp::upper ? inf : -inf );
  const AbstractLinAlgPack::SpVectorSlice::difference_type o = sv.offset();
  for( AbstractLinAlgPack::SpVectorSlice::const_iterator itr = sv.begin();
      itr != sv.end(); ++itr )
  {
    (*r)(itr->indice() + o) = itr->value();
  }
  DenseLinAlgPack::Vp_StV( r, -1.0, v );
  if( sign < 0 )
    DenseLinAlgPack::Vt_S( r, -1.0 );
}

#endif // 0
