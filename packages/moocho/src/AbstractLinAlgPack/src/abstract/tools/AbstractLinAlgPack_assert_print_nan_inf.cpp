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

#include <ostream>
#include <iomanip>

#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "RTOp_ROp_find_nan_inf.h"
#include "RTOpPack_RTOpC.hpp"
#include "check_nan_inf.h"
#include "Teuchos_TestForException.hpp"

namespace {

// Find a NaN or Inf element!
static RTOpPack::RTOpC                               find_nan_inf_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  find_nan_inf_targ;

class init_rtop_server_t {
public:
  init_rtop_server_t() {
    TEST_FOR_EXCEPT(0!=RTOp_ROp_find_nan_inf_construct(&find_nan_inf_op.op() ));
    find_nan_inf_targ = find_nan_inf_op.reduct_obj_create();
  }
}; 

init_rtop_server_t  init_rtop_server;

} // end namespace

bool AbstractLinAlgPack::assert_print_nan_inf( const value_type& val, char name[]
  , bool throw_excpt, std::ostream* out )
{
  if( RTOp_is_nan_inf(val) ) {
    std::ostringstream omsg;
    omsg
      << "The scalar \"" << name
      << "\" = " << val << " is not a valid bounded number";
    if(out)
      *out << omsg.str() << std::endl;
    TEST_FOR_EXCEPTION(
      throw_excpt,NaNInfException
      ,"assert_print_nan_inf(...) : Error, " << omsg.str() );
    return false;
  }
  return true;
}

bool AbstractLinAlgPack::assert_print_nan_inf(
  const Vector& v, char name[]
  ,bool throw_excpt, std::ostream* out
  )
{
  find_nan_inf_op.reduct_obj_reinit(&*find_nan_inf_targ);
  const Vector* vecs[1] = { &v };
  apply_op(find_nan_inf_op,1,vecs,0,NULL,&*find_nan_inf_targ);
  RTOp_ROp_find_nan_inf_reduct_obj_t
    ele =RTOp_ROp_find_nan_inf_val(find_nan_inf_op(*find_nan_inf_targ));
  if(out && ele.i) {
    *out
      << "The vector \"" << name << "\" has the first following NaN or Inf element\n"
      << name << "(" << ele.i << ") = " << ele.v0_i << std::endl;
  }
  TEST_FOR_EXCEPTION(
    ele.i && throw_excpt, NaNInfException
    ,"assert_print_nan_inf(...) : Error, the vector named "
    << name << " has at least one element which is NaN or Inf" );
  
  return ele.i == 0;
}
