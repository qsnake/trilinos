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

#include "MoochoPack_InitFinDiffReducedHessian_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 4;

  enum local_EOptions {
    INITIALIZATION_METHOD,
    MAX_COND,
    MIN_DIAG,
    STEP_SCALE
  };

  const char* local_SOptions[local_num_options]	= {
    "initialization_method",
    "max_cond",
    "min_diag",
    "step_scale"
  };

}

namespace MoochoPack {

InitFinDiffReducedHessian_StepSetOptions::InitFinDiffReducedHessian_StepSetOptions(
        InitFinDiffReducedHessian_Step* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      InitFinDiffReducedHessian_Step >( target )
{}

void InitFinDiffReducedHessian_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef InitFinDiffReducedHessian_Step target_t;
  switch( (local_EOptions)option_num ) {
      case INITIALIZATION_METHOD:
    {
      const std::string &option = option_value.c_str();
      if( option == "SCALE_IDENTITY" )
        target().initialization_method( target_t::SCALE_IDENTITY );
      else if( option == "SCALE_DIAGONAL" )
        target().initialization_method( target_t::SCALE_DIAGONAL );
      else if( option == "SCALE_DIAGONAL_ABS" )
        target().initialization_method( target_t::SCALE_DIAGONAL_ABS );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"initialization_method\"." );
      break;
    }
      case MAX_COND:
      target().max_cond(std::fabs(std::atof(option_value.c_str())));
      break;
    case MIN_DIAG:
      target().min_diag(std::abs(std::atoi(option_value.c_str())));
      break;
    case STEP_SCALE:
      target().step_scale(std::fabs(std::atof(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
