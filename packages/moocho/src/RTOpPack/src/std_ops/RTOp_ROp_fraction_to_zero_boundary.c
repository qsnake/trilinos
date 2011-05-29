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
/* RTOp_ROp_fraction_to_zero_boundary.c */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/8/2002 at 19:19 */
/* */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )


#include "RTOp_ROp_fraction_to_zero_boundary.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for operator object instance data */
#include "RTOp_reduct_min_value.h"  /* Reduction of intermediate reduction objects */



/* Implementation functions for RTOp_RTOp */

static int RTOp_ROp_fraction_to_zero_boundary_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  /* */
  /* Declare local variables */
  /* */

  /* Access to the operator object instance data */
  RTOp_value_type *tau = (RTOp_value_type*)obj_data;
  /* Access to the reduction object data */
  RTOp_value_type *alpha_max = (RTOp_value_type*)reduct_obj;
  /* Vector data */
  RTOp_index_type           sub_dim;
  /* v0 */
  const RTOp_value_type     *v0_val;
  ptrdiff_t                 v0_val_s;
  /* v1 */
  const RTOp_value_type     *v1_val;
  ptrdiff_t                 v1_val_s;

  /* Automatic temporary variables */
  register RTOp_index_type  k;
  /* Temporary element-wise reduction object */
  RTOp_value_type alpha_max_ith;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 2 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( /* Validate sub_dim */
    vecs[1].sub_dim != vecs[0].sub_dim
    )
    return RTOp_ERR_INCOMPATIBLE_VECS;
  assert(obj_data);
  assert(reduct_obj);

  /* */
  /* Get pointers to data */
  /* */
  sub_dim       = vecs[0].sub_dim;
  /* v0 */
  v0_val        = vecs[0].values;
  v0_val_s      = vecs[0].values_stride;
  /* v1 */
  v1_val        = vecs[1].values;
  v1_val_s      = vecs[1].values_stride;

  /* */
  /* Apply the operator: */
  /* */
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s, v1_val += v1_val_s )
  {
    /* Element-wise reduction */
    alpha_max_ith = ((*v1_val) >= 0) ? 1.0 : (*tau)*(*v0_val)/(-(*v1_val));
    /* Reduction of intermediates */
    (*alpha_max) = min( (*alpha_max), alpha_max_ith );
  }

  return 0; /* success? */
}

static int RTOp_ROp_fraction_to_zero_boundary_reduct_obj_reinit(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  ,RTOp_ReductTarget reduct_obj
  )
{
  RTOp_value_type *alpha_max = (RTOp_value_type*)reduct_obj;
  *alpha_max = 1.0;
  return 0;
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_fraction_to_zero_boundary_vtbl =
{
  &RTOp_obj_value_vtbl
  ,&RTOp_obj_value_vtbl
  ,"ROp_fraction_to_zero_boundary"
  ,RTOp_ROp_fraction_to_zero_boundary_reduct_obj_reinit
  ,RTOp_ROp_fraction_to_zero_boundary_apply_op
  ,RTOp_reduct_min_value
  ,RTOp_get_reduct_min_value_op
};

/* Class specific functions */

int RTOp_ROp_fraction_to_zero_boundary_construct( RTOp_value_type tau,  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_ROp_fraction_to_zero_boundary_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return RTOp_ROp_fraction_to_zero_boundary_init(tau,op);
}

int RTOp_ROp_fraction_to_zero_boundary_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}

int RTOp_ROp_fraction_to_zero_boundary_init( RTOp_value_type tau, struct RTOp_RTOp* op )
{
  RTOp_value_type *ptr_tau = (RTOp_value_type*)op->obj_data;
  *ptr_tau = tau;
  return 0;
}
RTOp_value_type RTOp_ROp_fraction_to_zero_boundary_val(RTOp_ReductTarget reduct_obj)
{
  return *((RTOp_value_type*)reduct_obj);
}

