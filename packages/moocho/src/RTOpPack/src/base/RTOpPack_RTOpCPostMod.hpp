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

// /////////////////////////////////////////////////////////////////////////
// RTOpPack_RTOpCPostMod.hpp

#ifndef RTOPPACK_RTOP_C_POST_MOD_HPP
#define RTOPPACK_RTOP_C_POST_MOD_HPP

#include "RTOpPack_RTOpC.hpp"

namespace RTOpPack {

class RTOpCPostMod {
public:

  /** \brief . */
  RTOpCPostMod( const RTOp_RTOp_vtbl_t *vtbl ) : vtbl_(vtbl)
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        !(vtbl && vtbl->obj_data_vtbl && vtbl->obj_data_vtbl->obj_create)
        ,std::logic_error, "Error!"
        );
#endif			
    }
  /** \brief . */
  void initialize(RTOpC *op) const
    {
      op->op().vtbl = vtbl_;
      op->op().vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->op().obj_data);
    }
  
private:
  
  const RTOp_RTOp_vtbl_t *vtbl_;
  
  RTOpCPostMod(); // Not defined and not to be called.

};

} // namespace RTOpPack

#endif // RTOPPACK_RTOP_C_POST_MOD_HPP
