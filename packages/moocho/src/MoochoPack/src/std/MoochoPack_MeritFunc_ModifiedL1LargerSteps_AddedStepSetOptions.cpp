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

#include "MoochoPack_MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 5;

  const char options_group_name[] = "MeritFuncModifiedL1LargerSteps";

  enum local_EOptions {
    AFTER_K_ITER,
    OBJ_INCREASE_THRESHOLD,
    MAX_POS_PENALTY_INCREASE,
    POS_TO_NEG_PENALTY_INCREASE,
    INCR_MULT_FACTOR
  };

  const char* local_SOptions[local_num_options]	= {
    "after_k_iter",
    "obj_increase_threshold",
    "max_pos_penalty_increase",
    "pos_to_neg_penalty_increase",
    "incr_mult_factor"
  };

}

namespace MoochoPack {

MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions::MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions(
      MeritFunc_ModifiedL1LargerSteps_AddedStep* target )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        options_group_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      MeritFunc_ModifiedL1LargerSteps_AddedStep >( target )
{}

void MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  switch( (local_EOptions)option_num ) {
    case AFTER_K_ITER: {
      target().after_k_iter( ::atoi( option_value.c_str() ) );
      break;
    }
    case OBJ_INCREASE_THRESHOLD: {
      target().obj_increase_threshold( ::atof( option_value.c_str() ) );
      break;
    }
    case MAX_POS_PENALTY_INCREASE: {
      target().max_pos_penalty_increase( ::fabs( ::atof( option_value.c_str() ) ) );
      break;
    }
    case POS_TO_NEG_PENALTY_INCREASE: {
      target().pos_to_neg_penalty_increase( ::fabs( ::atof( option_value.c_str() ) ) );
      break;
    }
    case INCR_MULT_FACTOR: {
      target().incr_mult_factor( ::fabs( ::atof( option_value.c_str() ) ) );
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 

#endif // 0
