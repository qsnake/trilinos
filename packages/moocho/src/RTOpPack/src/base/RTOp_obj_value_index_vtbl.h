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

#ifndef RTOP_OBJ_VALUE_INDEX_VTBL_H
#define RTOP_OBJ_VALUE_INDEX_VTBL_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Virtual function table for a simple object of a <tt>{RTOp_value_type,RTOp_index_type}</tt> pair.
  *
  * The functions do the following:
  * <ul>
  *	<li> <tt>get_obj_type_num_entries(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>num_values</tt>    [out] Returns 1
  *		<li> <tt>num_indexes</tt>   [out] Returns 1
  *		<li> <tt>num_chars</tt>     [out] Returns 0
  *		</ul>
  *	<li> <tt>obj_create(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [out] Points an allocated to a <tt>RTOp_value_index_type</tt>
  *                                 object initialized to <tt>{0.0,0}</tt>
  *		</ul>
  *	<li> <tt>obj_reinit(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [in/out] <tt>RTOp_value_index_type</tt> object reinitialized to <tt>{0.0,0}</tt>
  *		</ul>
  *	<li> <tt>obj_free(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [in/out] allocated object is freed and <tt>obj</tt>
  *                                 is set to NULL.
  *		</ul>
  *	<li> <tt>extract_state(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *     <li> <tt>obj</tt>           [in] Allocated <tt>RTOp_value_index_type</tt> object.
  *		<li> <tt>num_values</tt>    [in] Must be 1
  *     <li> <tt>value_data</tt>    [out] <tt>value_data[0] = ((RTOp_value_index_type*)obj->value)</tt>
  *		<li> <tt>num_indexes</tt>   [in] Must be 1
  *     <li> <tt>index_data</tt>    [out] <tt>index_data[0] = ((RTOp_value_index_type*)obj->index)</tt>
  *		<li> <tt>num_chars</tt>     [in] Must be 0
  *     <li> <tt>char_data</tt >    [out] Must be NULL
  *		</ul>
  *	<li> <tt>load_state(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>num_values</tt>    [in] Must be 1
  *     <li> <tt>value_data</tt>    [in] <tt>value_data[0]</tt>
  *		<li> <tt>num_indexes</tt>   [in] Must be 1
  *     <li> <tt>index_data</tt>    [in] <tt>index_data[0]</tt>
  *		<li> <tt>num_chars</tt>     [in] Must be 0
  *     <li> <tt>char_data</tt >    [in] Must be NULL
  *     <li> <tt>obj</tt>           [in/out] If <tt>*obj == NULL</tt> then
  *                                 <tt>*obj = malloc(sizeof(RTOp_value_index_type))</tt> and
  *                                 <tt>*(RTOp_value_index_type)*obj = {value_data[0],index_data[0]}</tt>
  *		</ul>
  * </ul>
  */
/*@{ */

/* Virtual function table */
extern const struct RTOp_obj_type_vtbl_t   RTOp_obj_value_index_vtbl;

/* Object type structure for a value, index pair. */
struct RTOp_value_index_type {
	/* */
	RTOp_value_type  value;
	/* */
	RTOp_index_type  index;
};

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_OBJ_VALUE_INDEX_VTBL_H */
