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

#ifndef RTOP_GET_REDUCT_OP_H
#define RTOP_GET_REDUCT_OP_H

#define INSERT_GET_REDUCT_OP_FUNCS(                                                                                   \
    num_values,num_indexes,num_chars,reduct_obj_t,reduce_reduct_obj,targ_load_state,targ_extract_state                \
    ,external_reduct_op,get_reduct_op                                                                                 \
    )                                                                                                                 \
static void CALL_API external_reduct_op( void* in_targ_array, void* inout_targ_array                                           \
  , int* len, RTOp_Datatype* datatype )                                                                             \
{                                                                                                                     \
  struct reduct_obj_t                                                                                               \
    in_obj, *in_obj_p = &in_obj, inout_obj, *inout_obj_p = &inout_obj;                                            \
    const int                                                                                                         \
        values_off  = 3*sizeof(RTOp_value_type),                                                                      \
        indexes_off = values_off + num_values*sizeof(RTOp_value_type),                                                \
        chars_off   = indexes_off + num_indexes*sizeof(RTOp_index_type);                                              \
  const int size_obj = chars_off + (num_chars)*sizeof(RTOp_index_type);                                             \
  char                                                                                                              \
     *in_array    = in_targ_array,                                                                                 \
    *inout_array = inout_targ_array;                                                                              \
  int i;                                                                                                            \
  for( i = 0; i < *len; ++i, in_array += size_obj, inout_array += size_obj ) {                                      \
    targ_load_state(                                                                                              \
      NULL, NULL                                                                                                \
      ,num_values,  num_values  ? (RTOp_value_type*)(in_array + values_off) : NULL                              \
      ,num_indexes, num_indexes ? (RTOp_index_type*)(in_array + values_off) : NULL                              \
      ,num_chars,   num_chars   ? (RTOp_char_type*) (in_array + values_off) : NULL                              \
      ,(void**)&in_obj_p );                                                                                     \
    targ_load_state(                                                                                              \
      NULL, NULL                                                                                                \
      ,num_values,  num_values  ? (RTOp_value_type*)(inout_array + values_off) : NULL                           \
      ,num_indexes, num_indexes ? (RTOp_index_type*)(inout_array + values_off) : NULL                           \
      ,num_chars,   num_chars   ? (RTOp_char_type*) (inout_array + values_off) : NULL                           \
      ,(void**)&inout_obj_p );                                                                                  \
    reduce_reduct_objs( NULL, NULL, &in_obj, &inout_obj );                                                        \
    targ_extract_state(                                                                                           \
      NULL, NULL, &inout_obj                                                                                    \
      ,num_values,  num_values  ? (RTOp_value_type*)(inout_array + values_off) : NULL                           \
      ,num_indexes, num_indexes ? (RTOp_index_type*)(inout_array + values_off) : NULL                           \
      ,num_chars,   num_chars   ? (RTOp_char_type*) (inout_array + values_off) : NULL                           \
      );                                                                                                        \
  }                                                                                                                 \
}                                                                                                                     \
static int get_reduct_op(                                                                                             \
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data                                                         \
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )                                                                 \
{                                                                                                                     \
  *reduct_op_func_ptr = external_reduct_op;                                                                         \
  return 0;                                                                                                         \
}

#endif
