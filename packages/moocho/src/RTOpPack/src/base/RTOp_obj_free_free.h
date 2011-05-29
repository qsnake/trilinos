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

#ifndef RTOP_OBJ_FREE_FREE_H
#define RTOP_OBJ_FREE_FREE_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* */
/** Definition of function to be used for </tt>obj_free<tt> in RTOp_RTOp_vtbl_t which
 * simply calls </tt>free(...)<tt>.
 *
 * @param  vtbl           [in] Totally ignored.
 * @param  instance_data  [in] Totally ignored.
 * @param  obj            [in/out]  If </tt>*obj != NULL<tt> on input then </tt>free(*obj)<tt> is called.
 *                        On output, </tt>*obj<tt> is set to NULL.
 *
 * @return Returns </tt>0<tt> for success.
 */
int RTOp_obj_free_free( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void** obj );

#ifdef __cplusplus
}
#endif

#endif /* RTOP_OBJ_FREE_FREE_H */
