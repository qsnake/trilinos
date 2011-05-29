/*
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
*/

/* ////////////////////////////////////////////////////////////// */
/* RTOp_apply_op_serial.c */

#include "RTOp_apply_op_serial.h"

#include <stdlib.h>

int RTOp_apply_op_serial(
  RTOp_index_type full_dim
  ,const int      num_vecs,  const RTOp_value_type*      vec_ptrs[],  const ptrdiff_t      vec_strides[]
  ,const int num_targ_vecs,  RTOp_value_type*       targ_vec_ptrs[],  const ptrdiff_t targ_vec_strides[]
  ,const RTOp_index_type first_ele_in, const RTOp_index_type sub_dim_in, const RTOp_index_type global_offset_in
  ,const struct RTOp_RTOp* op
  ,RTOp_ReductTarget reduct_obj
  )
{
  int                          err            = 0;
  RTOp_index_type              sub_dim        = 0;
  struct RTOp_SubVector        *sub_vecs      = NULL;
  struct RTOp_MutableSubVector *targ_sub_vecs = NULL;
  int                          k;
  /* Sort out the input and get the number of vector elements to operator over */
#ifdef RTOp_DEBUG
  assert( num_vecs || num_targ_vecs );
  if(num_vecs)
    assert( vec_ptrs != NULL );
  if(num_targ_vecs)
    assert( targ_vec_ptrs != NULL );
  assert( 0 <= sub_dim_in && sub_dim_in <= full_dim );
#endif
  sub_dim = sub_dim_in ? sub_dim_in : full_dim - (first_ele_in - 1); /* Dimension of logical vectors */
  /* Create the sub-vector data structures */
  if(num_vecs) {
    sub_vecs = malloc( sizeof(struct RTOp_SubVector) * num_vecs );
    for( k = 0; k < num_vecs; ++k ) {
#ifdef RTOp_DEBUG
      assert( vec_ptrs[k] != NULL );
#endif
      RTOp_sub_vector(
        global_offset_in
        ,sub_dim
        ,vec_ptrs[k] + (first_ele_in -1) * vec_strides[k]
        ,vec_strides[k]
        ,&sub_vecs[k]
        );
    }
  }
  if(num_targ_vecs) {
    targ_sub_vecs = malloc( sizeof(struct RTOp_MutableSubVector) * num_targ_vecs );
    for( k = 0; k < num_targ_vecs; ++k ) {
#ifdef RTOp_DEBUG
      assert( targ_vec_ptrs[k] != NULL );
#endif
      RTOp_mutable_sub_vector(
        global_offset_in
        ,sub_dim
        ,targ_vec_ptrs[k] + (first_ele_in -1) * targ_vec_strides[k]
        ,targ_vec_strides[k]
        ,&targ_sub_vecs[k]
        );
    }
  }
  /* Apply the reduction/transformation operator in one chunk */
  err = RTOp_apply_op( op, num_vecs, sub_vecs, num_targ_vecs, targ_sub_vecs, reduct_obj );
  /* Free the sub-vector data structures */
  if(      sub_vecs ) free(      sub_vecs );
  if( targ_sub_vecs ) free( targ_sub_vecs );

  return err;  /* This could be an error code! */
}
