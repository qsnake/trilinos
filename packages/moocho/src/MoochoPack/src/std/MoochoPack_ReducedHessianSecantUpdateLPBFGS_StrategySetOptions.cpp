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

#include "MoochoPack_ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 4;

  enum local_EOptions {
    MIN_NUM_UPDATES_PROJ_START
    ,MAX_NUM_UPDATES_PROJ_START
    ,NUM_SUPERBASICS_SWITCH_DENSE
    ,NUM_ADD_RECENT_UPDATES
  };

  const char* local_SOptions[local_num_options]	= {
    "min_num_updates_proj_start"
    ,"max_num_updates_proj_start"
    ,"num_superbasics_switch_dense"
    ,"num_add_recent_updates"
  };

}

namespace MoochoPack {

ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::ReducedHessianSecantUpdateLPBFGS_StrategySetOptions(
  ReducedHessianSecantUpdateLPBFGS_Strategy* target
  , const char opt_grp_name[] )
  : OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions )
  , OptionsFromStreamPack::SetOptionsToTargetBase< ReducedHessianSecantUpdateLPBFGS_Strategy >( target )
{}

void ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::setOption(
  int option_num, const std::string& option_value )
{
  switch( (local_EOptions)option_num ) {
    case MIN_NUM_UPDATES_PROJ_START: {
      target().min_num_updates_proj_start( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    case MAX_NUM_UPDATES_PROJ_START: {
      target().max_num_updates_proj_start( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    case NUM_SUPERBASICS_SWITCH_DENSE: {
      target().num_superbasics_switch_dense( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    case NUM_ADD_RECENT_UPDATES: {
      target().num_add_recent_updates( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 

#endif // 0
