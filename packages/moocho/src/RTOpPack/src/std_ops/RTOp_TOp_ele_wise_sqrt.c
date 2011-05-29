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

#include <math.h>

/* ///////////////////////////////////////////// */
/* RTOp_TOp_ele_wise_sqrt.c */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/2/2002 at 19:1 */
/* */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_TOp_ele_wise_sqrt.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for operator object instance data */



/* Implementation functions for RTOp_RTOp */

static int RTOp_TOp_ele_wise_sqrt_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  /* */
  /* Declare local variables */
  /* */

  /* Vector data */
  RTOp_index_type           sub_dim;
  /* z0 */
  RTOp_value_type           *z0_val;
  ptrdiff_t                 z0_val_s;

  /* Automatic temporary variables */
  register RTOp_index_type  k;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 0 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;


  /* */
  /* Get pointers to data */
  /* */
  sub_dim       = targ_vecs[0].sub_dim;
  /* z0 */
  z0_val        = targ_vecs[0].values;
  z0_val_s      = targ_vecs[0].values_stride;


  /* */
  /* Apply the operator: */
  /* */
  for( k = 0; k < sub_dim; ++k, z0_val += z0_val_s )
    {
    /* Element-wise transformation */
    assert(*z0_val >= 0);
    (*z0_val) = sqrt((*z0_val));
    }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_ele_wise_sqrt_vtbl =
{
  &RTOp_obj_null_vtbl
  ,&RTOp_obj_null_vtbl
  ,"TOp_ele_wise_sqrt"
  ,NULL
  ,RTOp_TOp_ele_wise_sqrt_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_ele_wise_sqrt_construct(  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_ele_wise_sqrt_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return 0;
}

int RTOp_TOp_ele_wise_sqrt_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}



