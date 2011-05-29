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

#ifndef RTOP_OBJ_VALUES_VTBL_H
#define RTOP_OBJ_VALUES_VTBL_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* */
/** Virtual function table for a simple array of scalar objects of type <tt>RTOp_value_type</tt>.
 *
 * The functions do the following:
 * <ul>
 *	<li> <tt>get_obj_type_num_entries(...)</tt>
 *		<ul>
 *		<li> <tt>vtbl</tt>          [in] Ignored
 *		<li> <tt>instance_data</tt> [in] Must be <tt>!=NULL</tt> and points to an \c RTOp_index_type
 *                                  that gives the number of values in the list.
 *		<li> <tt>num_values</tt>    [out] Returns \c *(RTOp_index_type*)instance_data
 *		<li> <tt>num_indexes</tt>   [out] Returns 0
 *		<li> <tt>num_chars</tt>     [out] Returns 0
 *		</ul>
 *	<li> <tt>obj_create(...)</tt>
 *		<ul>
 *		<li> <tt>vtbl</tt>          [in] Ignored
 *		<li> <tt>instance_data</tt> [in] Must be <tt>!=NULL</tt> and points to an \c RTOp_index_type
 *                                  that gives the number of values in the list.
 *		<li> <tt>obj</tt>           [out] Points an allocated list of <tt>RTOp_value_type</tt> values
 *                                  all initialized to 0.0.
 *		</ul>
 *	<li> <tt>obj_reinit(...)</tt>
 *		<ul>
 *		<li> <tt>vtbl</tt>          [in] Ignored
 *		<li> <tt>instance_data</tt> [in] Must be <tt>!=NULL</tt> and points to an \c RTOp_index_type
 *                                  that gives the number of values in the list.
 *		<li> <tt>obj</tt>           [in/out] Points an allocated list of <tt>RTOp_value_type</tt> values
 *                                  all initialized to 0.0.
 *		</ul>
 *	<li> <tt>obj_free(...)</tt>
 *		<ul>
 *		<li> <tt>vtbl</tt>          [in] Ignored
 *		<li> <tt>instance_data</tt> [in] Must be <tt>!=NULL</tt> and points to an \c RTOp_index_type
 *                                  that gives the number of values in the list.
 *		<li> <tt>obj</tt>           [in/out] allocated object is freed and <tt>obj</tt>
 *                                  is set to NULL.
 *		</ul>
 *	<li> <tt>extract_state(...)</tt>
 *		<ul>
 *		<li> <tt>vtbl</tt>          [in] Ignored
 *		<li> <tt>instance_data</tt> [in] Must be <tt>!=NULL</tt> and points to an \c RTOp_index_type
 *                                  that gives the number of values in the list.
 *      <li> <tt>obj</tt>           [in] Allocated <tt>RTOp_index_type</tt> object.
 *		<li> <tt>num_values</tt>    [in] Must be \c *(RTOp_index_type*)instance_data
 *      <li> <tt>value_data</tt>    [out] <tt>value_data[k] = ((RTOp_value_type*)obj)[k], k = 1 ... num_values</tt>
 *		<li> <tt>num_indexes</tt>   [in] Must be 0
 *      <li> <tt>index_data</tt>    [out] Must be NULL
 *		<li> <tt>num_chars</tt>     [in] Must be 0
 *      <li> <tt>char_data</tt >    [out] Must be NULL
 *		</ul>
 *	<li> <tt>load_state(...)</tt>
 *		<ul>
 *		<li> <tt>vtbl</tt>          [in] Ignored
 *		<li> <tt>instance_data</tt> [in] Must be <tt>!=NULL</tt> and points to an \c RTOp_index_type
 *                                  that gives the number of values in the list.
 *		<li> <tt>num_values</tt>    [in] Must be \c *(RTOp_index_type*)instance_data
 *		<li> <tt>num_indexes</tt>   [in] Must be 0
 *      <li> <tt>index_data</tt>    [in] Must be NULL
 *		<li> <tt>num_chars</tt>     [in] Must be 0
 *      <li> <tt>char_data</tt >    [in] Must be NULL
 *      <li> <tt>obj</tt>           [in/out] If <tt>*obj == NULL</tt> then the object is allocated as in
 *                                  \c obj_create() before it is assigned.  On output:
 *                                  <tt>((RTOp_value_type*)obj)[k] = value_data[k], k = 1 ... num_values</tt>
 *		</ul>
 *	</ul>
 */
extern const struct RTOp_obj_type_vtbl_t   RTOp_obj_values_vtbl;

#ifdef __cplusplus
}
#endif

#endif /* RTOP_OBJ_VALUES_VTBL_H */
