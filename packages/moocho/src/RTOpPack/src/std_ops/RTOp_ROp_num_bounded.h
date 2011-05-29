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

#ifndef RTOP_ROP_NUM_BOUNDED_H
#define RTOP_ROP_NUM_BOUNDED_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_num_bounded.h Reduction operator counts the number of elements with finite bounds.
  *
  * <tt>targ_obj <- size( { i | v[0](i) > -inf_bnd || v[1](i) < +inf_bnd } )</tt>
  *
  * This operator is defined to allow exactly two vector arguments
  * (<tt>num_vecs == 2</tt>) <tt>v[0]</tt>, <tt>v[1]</tt> and can only handle dense vectors.
  */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_num_bounded_vtbl;

/* Constructor */
int RTOp_ROp_num_bounded_construct( RTOp_value_type inf_bnd, struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_num_bounded_destroy( struct RTOp_RTOp* op );

/* Reset inf_bnd */
int RTOp_ROp_num_bounded_set_inf_bnd( RTOp_value_type inf_bnd, struct RTOp_RTOp* op );

/* Extract the number of bounded variables */
RTOp_index_type RTOp_ROp_num_bounded_val(RTOp_ReductTarget targ_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_NUM_BOUNDED_H */
