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

/* ///////////////////////////////////////////// */
/* RTOp_TOp_sign.h */

#ifndef RTOp_TOp_sign_H
#define RTOp_TOp_sign_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_sign.h
 *
 * This operator computes the sign of each element in v0 as:
 \verbatim
  z0(i) = sign(v0(i)), for i = 1...n
 \endverbatim
 * where
 \verbatim
                /  -1.0 : if alpha  < 0.0
 sign(alpha) =  |   0.0 : if alpha == 0.0
                \  +1.0 : if alpha  > 0.0
 \endverbatim
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_sign_vtbl;

/* Constructor */
int RTOp_TOp_sign_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_sign_destroy( struct RTOp_RTOp* op );




/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_sign_H */
