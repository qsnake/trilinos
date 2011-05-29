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

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 6/24/2002 at 21:2 */
/* */

#ifndef RTOp_ROp_max_rel_step_H
#define RTOp_ROp_max_rel_step_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_max_rel_step.h
 *
 * This reduction operator finds the maximum relative step change in \c v0 for:
 \verbatim

 v0 = v0 + v1;
 \endverbatim
 *
 * This is found by the reduciton:
 \verbatim

 element-wise reduction      : gamma = max{ |v1(i)| / ( 1.0 + |v0| ), i = 1...n }
 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_rel_step_vtbl;

/* Constructor */
int RTOp_ROp_max_rel_step_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_max_rel_step_destroy( struct RTOp_RTOp* op );

/* Extract the value of the reduction object \c gamma */
RTOp_value_type RTOp_ROp_max_rel_step_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_max_rel_step_H */
