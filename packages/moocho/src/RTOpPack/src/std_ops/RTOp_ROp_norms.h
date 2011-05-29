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

#ifndef RTOP_ROP_NORMS_H
#define RTOP_ROP_NORMS_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_norms.h Reduction operator classes for common norms.
  */

/** @name One norm reduction operator class.
 *
 * <tt>||v[0]||_1 -> targ_obj</tt>
 */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_1_vtbl;

/* Constructor */
int RTOp_ROp_norm_1_construct( struct RTOp_RTOp* op );

/* Extract the value of the norm */
RTOp_value_type RTOp_ROp_norm_1_val(RTOp_ReductTarget targ_obj);

/** @name Two (Euclidean) norm reduction operator class.
 *
 * <tt>||v[0]||_2 -> targ_obj</tt>
 */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_2_vtbl;

/* Constructor */
int RTOp_ROp_norm_2_construct( struct RTOp_RTOp* op );

/* Extract the value of the norm */
RTOp_value_type RTOp_ROp_norm_2_val(RTOp_ReductTarget targ_obj);

/** @name Infinity norm reduction operator class.
 *
 * <tt>||v[0]||_inf -> targ_obj</tt>
 */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_inf_vtbl;

/* Constructor */
int RTOp_ROp_norm_inf_construct( struct RTOp_RTOp* op );

/* Extract the value of the norm */
RTOp_value_type RTOp_ROp_norm_inf_val(RTOp_ReductTarget targ_obj);

/* Destructor (for all three norms) */
int RTOp_ROp_norm_destroy( struct RTOp_RTOp* op );

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_NORMS_H */
