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

#include "RTOp_apply_op_mpi.h"
#include "RTOp_parallel_helpers.h"
#include "RTOpToMPI.h"

#include <stdlib.h>

int RTOp_apply_op_mpi(
  MPI_Comm comm
  ,RTOp_index_type global_dim_in, RTOp_index_type local_sub_dim_in, RTOp_index_type local_offset_in
  ,const int num_cols
  ,const int      num_vecs,  const RTOp_value_type*      l_vec_ptrs[],  const ptrdiff_t      l_vec_strides[], const ptrdiff_t      l_vec_leading_dim[]
  ,const int num_targ_vecs,  RTOp_value_type*       l_targ_vec_ptrs[],  const ptrdiff_t l_targ_vec_strides[], const ptrdiff_t l_targ_vec_leading_dim[]
  ,const RTOp_index_type first_ele_in, const RTOp_index_type sub_dim_in, const RTOp_index_type global_offset_in
  ,const struct RTOp_RTOp* op
  ,RTOp_ReductTarget reduct_objs[]
  )
{
  int                          err                      = 0;
  struct RTOp_SubVector        *local_vecs              = NULL;
  struct RTOp_MutableSubVector *local_targ_vecs         = NULL;
  RTOp_index_type              overlap_first_local_ele  = 0;
  RTOp_index_type              overalap_local_sub_dim   = 0;
  RTOp_index_type              overlap_global_offset    = 0;
  int                          k;
  int                          kc;
  /* Validate the input */
#ifdef RTOp_DEBUG
  assert( num_vecs || num_targ_vecs );
  if(num_vecs)
    assert( l_vec_ptrs != NULL );
  if(num_targ_vecs)
    assert( l_targ_vec_ptrs != NULL );
  assert( 0 <= sub_dim_in && sub_dim_in <= global_dim_in );
#endif
  /* Pre-initialize the local sub-vectors */
  if(num_vecs) {
    local_vecs = malloc( sizeof(struct RTOp_SubVector) * num_vecs * num_cols );
    for( kc = 0; kc < num_cols; ++kc ) {
      for( k = 0; k < num_vecs; ++k )
        RTOp_sub_vector_null(&local_vecs[kc*num_cols+k]);
    }
  }
  if(num_targ_vecs) {
    local_targ_vecs = malloc( sizeof(struct RTOp_MutableSubVector) * num_targ_vecs );
    for( kc = 0; kc < num_cols; ++kc ) {
      for( k = 0; k < num_targ_vecs; ++k )
        RTOp_mutable_sub_vector_null(&local_targ_vecs[kc*num_cols+k]);
    }
  }
  /* Get the overlap in the current process with the input logical sub-vector */
  /* from (first_ele_in,sub_dim_in,global_offset_in) */
  RTOp_parallel_calc_overlap(
    global_dim_in, local_sub_dim_in, local_offset_in, first_ele_in, sub_dim_in, global_offset_in
    ,&overlap_first_local_ele, &overalap_local_sub_dim, &overlap_global_offset
    );
  if( overlap_first_local_ele != 0 ) {
    /* Sub-vector structs for the local elements that are to participate in the */
    /* reduction/transforamtion operation. */
    for( kc = 0; kc < num_cols; ++kc ) {
      for(k = 0; k < num_vecs; ++k) {
        RTOp_sub_vector(
          overlap_global_offset                                                  /* global_offset */
          ,overalap_local_sub_dim                                                /* sub_dim */
          ,l_vec_ptrs[k]+(overlap_first_local_ele-1)*l_vec_strides[k]
          + ( num_cols > 1 ? kc*l_vec_leading_dim[k] : 0 )                       /* values */
          ,l_vec_strides[k]                                                      /* values_stride */
          ,&local_vecs[kc*num_cols+k]
          );
      }
      for(k = 0; k < num_targ_vecs; ++k) {
        RTOp_mutable_sub_vector(
          overlap_global_offset                                                  /* global_offset */
          ,overalap_local_sub_dim                                                /* sub_dim */
          ,l_targ_vec_ptrs[k]+(overlap_first_local_ele-1)*l_targ_vec_strides[k]
          + ( num_cols > 1 ? kc*l_targ_vec_leading_dim[k] : 0 )                  /* values */
          ,l_targ_vec_strides[k]                                                 /* values_stride */
          ,&local_targ_vecs[kc*num_cols+k]
          );
      }
    }
  }
  /* */
  /* Apply the reduction operation over the sub-vectors in */
  /* this process then collect the reductions over */
  /* all the processes and return the result */
  /* to all the processes (including this one of course). */
  /* If all of the sub-svectors are empty then this will */
  /* just call the reduction operation with NULL sub-vectors */
  /* */
  err = RTOp_MPI_apply_op(
    comm, op, -1 /* MPI_Allreduce(...) */
    ,num_cols
    ,num_vecs,      num_vecs && overlap_first_local_ele      ? &local_vecs[0]      : NULL
    ,num_targ_vecs, num_targ_vecs && overlap_first_local_ele ? &local_targ_vecs[0] : NULL
    ,reduct_objs
    );

  if(local_vecs)      free(local_vecs);
  if(local_targ_vecs) free(local_targ_vecs);

  /* Deallocate memory */

  return err;
}
