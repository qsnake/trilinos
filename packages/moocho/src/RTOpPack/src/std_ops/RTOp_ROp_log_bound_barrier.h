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
/*       on 6/26/2002 at 21:9 */
/* */

#ifndef RTOp_ROp_log_bound_barrier_H
#define RTOp_ROp_log_bound_barrier_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_log_bound_barrier.h
 *
 * This operator computes the log barrier term for the doubly bounded
 * inequalities:
 \verbatim

 xl <= x < xu
 \endverbatim
 * The barrier is computed as:
 \verbatim

 sum{ log( x(i) - xl(i) ) + log( xu(i) - x(i) ) , for i = 1...n }
 \endverbatim
 * To call this operator you must pass in the vectors in the order:
 * <tt>v0 = x, v1 = xl, v2 = xu</tt>.
 *
 * The operation performed is:
 \verbatim

 element-wise reduction : log_result += log(v0(i) - v1(i)) + log(v2(i) - v0(i)), i = 1...n
 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_log_bound_barrier_vtbl;

/* Constructor */
int RTOp_ROp_log_bound_barrier_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_log_bound_barrier_destroy( struct RTOp_RTOp* op );


/* Extract the value of the reduction object <tt>log_result</tt> */
RTOp_value_type RTOp_ROp_log_bound_barrier_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_log_bound_barrier_H */
