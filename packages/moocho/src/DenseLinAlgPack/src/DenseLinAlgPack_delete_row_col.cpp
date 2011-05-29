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

#include <assert.h>

#include "DenseLinAlgPack_delete_row_col.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"

void DenseLinAlgPack::delete_row_col( size_type kd, DMatrixSliceTriEle* tri_M )
{
  // Validate input
  TEST_FOR_EXCEPT( !(  tri_M  ) );
  TEST_FOR_EXCEPT( !(  tri_M->rows()  ) );
  TEST_FOR_EXCEPT( !(  1 <= kd && kd <= tri_M->rows()  ) );

  DMatrixSlice   M = tri_M->gms();
  const size_type  n = M.rows();

  if( tri_M->uplo() == BLAS_Cpp::lower ) {
    // Move M31 up one row at a time
    if( 1 < kd && kd < n ) {
      Range1D rng(1,kd-1);
      for( size_type i = kd; i < n; ++i )
        M.row(i)(rng) = M.row(i+1)(rng);
    }
    // Move M33 up and to the left one column at a time
    if( kd < n ) {
      for( size_type i = kd; i < n; ++i )
        M.col(i)(i,n-1) = M.col(i+1)(i+1,n);
    }
  }
  else if(  tri_M->uplo() == BLAS_Cpp::upper ) {
    // Move M13 left one column at a time.
    if( 1 < kd && kd < n ) {
      Range1D rng(1,kd-1);
      for( size_type j = kd; j < n; ++j )
        M.col(j)(rng) = M.col(j+1)(rng);
    }
    // Move the updated U33 up and left one column at a time.
    if( kd < n ) {
      for( size_type j = kd; j < n; ++j )
        M.col(j)(kd,j) = M.col(j+1)(kd+1,j+1);
    }
  }
  else {
    TEST_FOR_EXCEPT(true); // Invalid input
  }
}
