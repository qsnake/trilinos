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

#include "MoochoPack_EvalNewPointTailoredApproach_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 1;

  enum local_EOptions {
    FD_DERIV_TESTING
  };

  const char* local_SOptions[local_num_options]	= {
    "fd_deriv_testing"
  };

}

namespace MoochoPack {

EvalNewPointTailoredApproach_StepSetOptions::EvalNewPointTailoredApproach_StepSetOptions(
        EvalNewPointTailoredApproach_Step* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      EvalNewPointTailoredApproach_Step >( target )
{}

void EvalNewPointTailoredApproach_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  typedef EvalNewPointTailoredApproach_Step target_t;
  switch( (local_EOptions)option_num ) {
      case FD_DERIV_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_DEFAULT" )
        target().fd_deriv_testing( target_t::FD_DEFAULT );
      else if( option == "FD_TEST" )
        target().fd_deriv_testing( target_t::FD_TEST );
      else if( option == "FD_NO_TEST" )
        target().fd_deriv_testing( target_t::FD_NO_TEST );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"fd_deriv_testing\".  Only the options "
          "FD_DEFAULT, FD_TEST, and FD_NO_TEST "
          "are available" );
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
