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

#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_obj_null_vtbl.h"

/* Implementation functions for RTOp_RTOp */

static int RTOp_TOp_assign_vectors_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_index_type        z_sub_dim;
  RTOp_value_type        *z_val = NULL;
  ptrdiff_t              z_val_s;
  RTOp_index_type        v0_sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;

  register RTOp_index_type k;
  RTOp_value_type *z_val_tmp = NULL;
#ifdef RTOp_DEBUG
  RTOp_index_type i;
#endif

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 1 || vecs == NULL )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || targ_vecs == NULL )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( targ_vecs[0].sub_dim != vecs[0].sub_dim )
    return RTOp_ERR_INCOMPATIBLE_VECS;

  /* */
  /* Get pointers to data */
  /* */

  /* z */
  z_sub_dim     = targ_vecs[0].sub_dim;
  z_val         = targ_vecs[0].values;
  z_val_s       = targ_vecs[0].values_stride;

  /* v0 */
  v0_sub_dim     = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  z_val_tmp = z_val;

  /* */
  /* Assign the elements */
  /* */

  if( z_val_s == 1 && v0_val_s == 1 ) {
    /* Slightly faster loop for unit stride vectors */
    for( k = 0; k < z_sub_dim; ++k )
      *z_val++ = *v0_val++;
  }
  else {
    /* More general implementation for one or both non-unit strides */
    for( k = 0; k < z_sub_dim; ++k, z_val += z_val_s, v0_val += v0_val_s )
      *z_val = *v0_val;
  }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_assign_vectors_vtbl =
{
  &RTOp_obj_null_vtbl  /* Use a null object for instance data */
  ,&RTOp_obj_null_vtbl /* use null type for target object */
  ,"TOp_assign_vectors"
  ,NULL /* use default from reduct_vtbl */
  ,RTOp_TOp_assign_vectors_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_assign_vectors_construct( struct RTOp_RTOp* op )
{
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_assign_vectors_vtbl;
  return 0; /* success? */
}

int RTOp_TOp_assign_vectors_destroy( struct RTOp_RTOp* op )
{
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0; /* success? */
}
