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

/* ///////////////////////////////////////////// */
/* RTOp_ROp_max.c */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/15/2002 at 17:33 */
/* */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_ROp_max.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for operator object instance data */
#include "RTOp_reduct_max_value.h"


/* Implementation functions for RTOp_RTOp */

static int RTOp_ROp_max_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  /* */
  /* Declare local variables */
  /* */

  /* Access to the reduction object data */
  RTOp_value_type *max_ele = (RTOp_value_type*)reduct_obj;
  /* Vector data */
  RTOp_index_type           sub_dim;
  /* v0 */
  const RTOp_value_type     *v0_val;
  ptrdiff_t                 v0_val_s;

  /* Automatic temporary variables */
  register RTOp_index_type  k;
  /* Temporary element-wise reduction object */
  RTOp_value_type max_ele_ith;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 1 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(reduct_obj);


  /* */
  /* Get pointers to data */
  /* */
  sub_dim       = vecs[0].sub_dim;
  /* v0 */
  v0_val        = vecs[0].values;
  v0_val_s      = vecs[0].values_stride;


  /* */
  /* Apply the operator: */
  /* */
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s )
  {
    /* Element-wise reduction */
    max_ele_ith = (*v0_val);
    /* Reduction of intermediates */
    (*max_ele) = max( (*max_ele), max_ele_ith );
  }

  return 0; /* success? */
}

static int RTOp_ROp_max_reduct_obj_reinit(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  ,RTOp_ReductTarget reduct_obj
  )
{
  RTOp_value_type *max_ele = (RTOp_value_type*)reduct_obj;
  *max_ele = -1e+50;
  return 0;
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_vtbl =
{
  &RTOp_obj_null_vtbl
  ,&RTOp_obj_value_vtbl
  ,"ROp_max"
  ,RTOp_ROp_max_reduct_obj_reinit
  ,RTOp_ROp_max_apply_op
  ,RTOp_reduct_max_value
  ,RTOp_get_reduct_max_value_op
};

/* Class specific functions */

int RTOp_ROp_max_construct(  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_ROp_max_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return 0;
}

int RTOp_ROp_max_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}


RTOp_value_type RTOp_ROp_max_val(RTOp_ReductTarget reduct_obj)
{
  return *((RTOp_value_type*)reduct_obj);
}

