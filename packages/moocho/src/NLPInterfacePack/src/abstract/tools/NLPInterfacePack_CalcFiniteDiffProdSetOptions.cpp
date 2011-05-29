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

#include <assert.h>
#include <math.h>

#include "NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

  const int local_num_options = 6;

  enum local_EOptions {
    FD_METHOD_ORDER
    ,FD_STEP_SELECT
    ,FD_STEP_SIZE
    ,FD_STEP_SIZE_MIN
    ,FD_STEP_SIZE_F
    ,FD_STEP_SIZE_C
  };

  const char* local_SOptions[local_num_options]	= {
    "fd_method_order"
    ,"fd_step_select"
    ,"fd_step_size"
    ,"fd_step_size_min"
    ,"fd_step_size_f"
    ,"fd_step_size_c"
  };

}

namespace NLPInterfacePack {

CalcFiniteDiffProdSetOptions::CalcFiniteDiffProdSetOptions(
  CalcFiniteDiffProd* target
  ,const char opt_grp_name[]
  )
  :OptionsFromStreamPack::SetOptionsFromStreamNode(opt_grp_name,local_num_options,local_SOptions)
  ,OptionsFromStreamPack::SetOptionsToTargetBase<CalcFiniteDiffProd>(target)
{}

void CalcFiniteDiffProdSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef CalcFiniteDiffProd target_t;
  switch( (local_EOptions)option_num ) {
      case FD_METHOD_ORDER:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_ORDER_ONE" )
        target().fd_method_order( target_t::FD_ORDER_ONE );
      else if( option == "FD_ORDER_TWO" )
        target().fd_method_order( target_t::FD_ORDER_TWO );
      else if( option == "FD_ORDER_TWO_CENTRAL" )
        target().fd_method_order( target_t::FD_ORDER_TWO_CENTRAL );
      else if( option == "FD_ORDER_TWO_AUTO" )
        target().fd_method_order( target_t::FD_ORDER_TWO_AUTO );
      else if( option == "FD_ORDER_FOUR" )
        target().fd_method_order( target_t::FD_ORDER_FOUR );
      else if( option == "FD_ORDER_FOUR_CENTRAL" )
        target().fd_method_order( target_t::FD_ORDER_FOUR_CENTRAL );
      else if( option == "FD_ORDER_FOUR_AUTO" )
        target().fd_method_order( target_t::FD_ORDER_FOUR_AUTO );
      else
        TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"CalcFiniteDiffProdSetOptions::setOption(...) : Error, incorrect value for "
          "\"fd_method_order\".  Only the options FD_ORDER_ONE, FD_ORDER_TWO, "
          "FD_ORDER_TWO_CENTRAL, FD_ORDER_TWO_AUTO, FD_ORDER_FOUR, FD_ORDER_FOUR_CENTRAL "
          "and FD_ORDER_FOUR_AUTO are available" );
      break;
    }
      case FD_STEP_SELECT:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_STEP_ABSOLUTE" )
        target().fd_step_select( target_t::FD_STEP_ABSOLUTE );
      else if( option == "FD_STEP_RELATIVE" )
        target().fd_step_select( target_t::FD_STEP_RELATIVE );
      else
        TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"CalcFiniteDiffProdSetOptions::setOption(...) : Error, incorrect value for "
          "\"fd_step_select\".  Only the options are available" );
      break;
    }
      case FD_STEP_SIZE:
      target().fd_step_size(std::atof(option_value.c_str()));
      break;
      case FD_STEP_SIZE_MIN:
      target().fd_step_size_min(std::atof(option_value.c_str()));
      break;
      case FD_STEP_SIZE_F:
      target().fd_step_size_f(std::atof(option_value.c_str()));
      break;
      case FD_STEP_SIZE_C:
      target().fd_step_size_c(std::atof(option_value.c_str()));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace NLPInterfacePack
