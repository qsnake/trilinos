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

#include "MoochoPack_TangentialStepWithInequStd_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

  const int local_num_options = 4;

  enum local_EOptions {
    WARM_START_FRAC
    ,QP_TESTING
    ,PRIMAL_FEASIBLE_POINT_ERROR
    ,DUAL_FEASIBLE_POINT_ERROR
  };

  const char* local_SOptions[local_num_options]	= {
    "warm_start_frac"
    ,"qp_testing"
    ,"primal_feasible_point_error"
    ,"dual_feasible_point_error"
  };

}

namespace MoochoPack {

TangentialStepWithInequStd_StepSetOptions::TangentialStepWithInequStd_StepSetOptions(
  TangentialStepWithInequStd_Step* target
  ,const char opt_grp_name[]
  )
  :OptionsFromStreamPack::SetOptionsFromStreamNode(
     opt_grp_name, local_num_options, local_SOptions )
  ,OptionsFromStreamPack::SetOptionsToTargetBase<
    TangentialStepWithInequStd_Step >( target )
{}

void TangentialStepWithInequStd_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  typedef TangentialStepWithInequStd_Step target_t;
  switch( (local_EOptions)option_num ) {
    case WARM_START_FRAC:
      target().warm_start_frac(std::fabs(std::atof(option_value.c_str())));
      break;
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
        TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"Error, incorrect value \'" << option << "\' for "
          "\"qp_testing\".  Only the options "
          "QP_TEST_DEFAULT, QP_TEST, and QP_NO_TEST "
          "are available" );
      break;
    }
        case PRIMAL_FEASIBLE_POINT_ERROR:
      target().primal_feasible_point_error(StringToBool("primal_feasible_point_error",option_value.c_str()));
      break;
      case DUAL_FEASIBLE_POINT_ERROR:
        target().dual_feasible_point_error(StringToBool("dual_feasible_point_error",option_value.c_str()));
      break;
      default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
