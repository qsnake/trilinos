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

#include "MoochoPack_FeasibilityStepReducedStd_StrategySetOptions.hpp"

// Define the options
namespace {

  const int local_num_options = 2;

  enum local_EOptions {
    QP_OBJECTIVE
    ,QP_TESTING
  };

  const char* local_SOptions[local_num_options]	= {
    "qp_objective"
    ,"qp_testing"
  };

}

namespace MoochoPack {

FeasibilityStepReducedStd_StrategySetOptions::FeasibilityStepReducedStd_StrategySetOptions(
        FeasibilityStepReducedStd_Strategy* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      FeasibilityStepReducedStd_Strategy >( target )
{}

void FeasibilityStepReducedStd_StrategySetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef FeasibilityStepReducedStd_Strategy target_t;
  switch( (local_EOptions)option_num ) {
      case QP_OBJECTIVE:
    {
      const std::string &option = option_value.c_str();
      if( option == "OBJ_MIN_FULL_STEP" )
        target().qp_objective( target_t::OBJ_MIN_FULL_STEP );
      else if( option == "OBJ_MIN_NULL_SPACE_STEP" )
        target().qp_objective( target_t::OBJ_MIN_NULL_SPACE_STEP );
      else if( option == "OBJ_RSQP" )
        target().qp_objective( target_t::OBJ_RSQP );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"qp_objective\".  Only the options "
          "OBJ_MIN_FULL_STEP, OBJ_MIN_NULL_SPACE_STEP, and  OBJ_RSQP"
          " are available" );
      break;
    }
      case QP_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "QP_TEST_DEFAULT" )
        target().qp_testing( target_t::QP_TEST_DEFAULT );
      else if( option == "QP_TEST" )
        target().qp_testing( target_t::QP_TEST );
      else if( option == "QP_NO_TEST" )
        target().qp_testing( target_t::QP_NO_TEST );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"qp_testing\".  Only the options "
          "QP_TEST_DEFAULT, QP_TEST, and QP_NO_TEST "
          "are available" );
      break;
    }
      default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
