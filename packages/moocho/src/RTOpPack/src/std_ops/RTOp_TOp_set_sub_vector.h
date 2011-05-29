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

#ifndef RTOP_TOP_SET_SUB_VECTOR_H
#define RTOP_TOP_SET_SUB_VECTOR_H

#include "RTOp.h"
#include "RTOp_SparseSubVector.h"

#ifdef __cplusplus
extern "C" {
#endif


/** \file RTOp_TOp_set_sub_vector.h Transforamtion operator for setting a sub-vector in a whole vector.
 *
 * <tt>z[0](sub_vec.global_offset+1,sub_vec.global_offset+sub_vec.sub_dim) <- sub_vec</tt>
 *
 * This operator is only defined to allow one vector argument
 * (<tt>num_targ_vecs == 1</tt>) <tt>z[0]</tt>.
 * Using a reduction operator to set a sub-vector will be reasonably
 * efficient for some types of vector subclasses (e.g. dense and sparse serial vectors
 * and parallel vectors with client in every process) but very slow for others
 * (e.g. out-of-core vectors and parallel vectors where operator object must be scattered
 * to various processes).
 * It would be better for vector subclasses to implement this operation directly
 * but if they don't you can use this operator to set the required sub-vector.
 *
 * This operator class works by allocating an internal sub-vector as the operator object
 * state data.
 *
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_set_sub_vector_vtbl;

/* */
/** Constructor.
 *
 * Note that a copy of sub_vec is not made.  Therefore, the client must
 * not disturb <tt>sub_vec</tt> while this operator object is in use!
 */
int RTOp_TOp_set_sub_vector_construct(
  const struct RTOp_SparseSubVector* sub_vec, struct RTOp_RTOp* op );

/* */
/** Reinitialize the range for the sub-vector to extract.
 *
 * Note that a copy of sub_vec is not made.  Therefore, the client must
 * not disturb <tt>sub_vec</tt> while this operator object is in use!
 */
int RTOp_TOp_set_sub_vector_set_sub_vec(
  const struct RTOp_SparseSubVector* sub_vec, struct RTOp_RTOp* op );

/* */
/** Destructor.
 */
int RTOp_TOp_set_sub_vector_destroy( struct RTOp_RTOp* op );

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_TOP_SET_SUB_VECTOR_H */
