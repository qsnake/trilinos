#if 0

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

#include "ReducedHessianBFGSStd_StepSetOptions.h"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 5;

  enum local_EOptions {
    RESCALE_INIT_IDENTITY
    ,USE_DAMPENING
    ,SECANT_TESTING
    ,SECANT_WARNING_TOL
    ,SECANT_ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
    "rescale_init_identity"
      ,"use_dampening"
    ,"secant_testing"
    ,"secant_warning_tol"
      ,"secant_error_tol"
  };

}

namespace MoochoPack {

ReducedHessianBFGSStd_StepSetOptions::ReducedHessianBFGSStd_StepSetOptions(
        ReducedHessianBFGSStd_Step* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      ReducedHessianBFGSStd_Step >( target )
{}

void ReducedHessianBFGSStd_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  typedef ReducedHessianBFGSStd_Step target_t;
  switch( (local_EOptions)option_num ) {
      case RESCALE_INIT_IDENTITY:
      target().rescale_init_identity(
        StringToBool( "rescale_init_identity", option_value.c_str() ));
      break;
      case USE_DAMPENING:
      target().use_dampening(
        StringToBool( "use_dampening", option_value.c_str() ));
      break;
      case SECANT_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "DEFAULT" )
        target().secant_testing( target_t::SECANT_TEST_DEFAULT );
      else if( option == "TEST" )
        target().secant_testing( target_t::SECANT_TEST_ALWAYS );
      else if( option == "NO_TEST" )
        target().secant_testing( target_t::SECANT_NO_TEST );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"secant_testing\"." );
      break;
    }
      case SECANT_WARNING_TOL:
      target().secant_warning_tol(::fabs(::atof(option_value.c_str())));
      break;
      case SECANT_ERROR_TOL:
      target().secant_error_tol(::fabs(::atof(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 

#endif // 0
