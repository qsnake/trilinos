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

#ifndef RTOP_TOP_SET_ELE_H
#define RTOP_TOP_SET_ELE_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_set_ele.h Set an individual element of a vector.
  *
  * <tt>targ_vec(i) <- alpha</tt>
  *
  * This transformation operator only sets an individual element of
  * a vector and leaves all of the others the same.
  *
  * This operator is only defined for a self transformation (<tt>num_vecs == 0</tt>).
  */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_set_ele_vtbl;

/* Constructor */
int RTOp_TOp_set_ele_construct( RTOp_index_type i, RTOp_value_type alpha
  , struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_set_ele_destroy( struct RTOp_RTOp* op );

/* Reset i and alpha */
int RTOp_TOp_set_ele_set_i_alpha( RTOp_index_type i, RTOp_value_type alpha
  , struct RTOp_RTOp* op );

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_TOP_SET_ELE_H */
