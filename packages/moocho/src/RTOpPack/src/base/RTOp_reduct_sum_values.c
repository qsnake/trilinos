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

#include "RTOp_reduct_sum_values.h"

int RTOp_reduct_sum_values(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget in_targ_obj, RTOp_ReductTarget inout_targ_obj )
{
  int num_values, k;
#ifdef RTOp_DEBUG
  assert(obj_data);
#endif
  num_values = *(RTOp_index_type*)obj_data;
  /* inout_dot_prod += in_dot_prod */
  for( k = 0; k < num_values; ++k )
    ((RTOp_value_type*)inout_targ_obj)[k] += ((RTOp_value_type*)in_targ_obj)[k];
  return 0;
}

static void CALL_API external_reduct_op( void* in_targ_array, void* inout_targ_array
  , int* len, RTOp_Datatype* datatype )
{
  /* inout_dot_prod += in_dot_prod */
  RTOp_index_type
    num_values   = *(RTOp_value_type*)in_targ_array; /* num_values is first size member */
  RTOp_value_type /* index past the size members */
    *in_targs    = (RTOp_value_type*)in_targ_array    + 3,
    *inout_targs = (RTOp_value_type*)inout_targ_array + 3;
  int i, k;
  for( i = 0; i < *len; ++i, inout_targs += (3 + num_values), in_targs += (3 + num_values) ) {
    for( k = 0; k < num_values; ++k )
      inout_targs[k] += in_targs[k];
  }
}

int RTOp_get_reduct_sum_values_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )
{
  *reduct_op_func_ptr = external_reduct_op;
  return 0;
}
