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

#include "MoochoPack_LineSearchFilter_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

const int local_num_options = 11;

enum local_EOptions 
  {
    GAMMA_THETA
    ,GAMMA_F
    ,F_MIN
    ,GAMMA_ALPHA
    ,DELTA
    ,S_THETA
    ,S_F
    ,THETA_SMALL_FACT
    ,THETA_MAX
    ,ETA_F
    ,BACK_TRACK_FRAC
  };

const char* local_SOptions[local_num_options] = 
  {
    "gamma_theta"
    ,"gamma_f"
    ,"f_min"
    ,"gamma_alpha"
    ,"delta"
    ,"s_theta"
    ,"s_f"
    ,"theta_small_fact"
    ,"theta_max"
    ,"eta_f"
    ,"back_track_frac"
  };

}

namespace MoochoPack {

LineSearchFilter_StepSetOptions::LineSearchFilter_StepSetOptions(
  LineSearchFilter_Step* target
  , const char opt_grp_name[] )
  :
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase< LineSearchFilter_Step >( target )
  {
  }

void LineSearchFilter_StepSetOptions::setOption( 
  int option_num, const std::string& option_value )
  {
  using OptionsFromStreamPack::StringToBool;
  
  typedef LineSearchFilter_Step target_t;
  switch( (local_EOptions)option_num ) {
    case GAMMA_THETA:
      target().gamma_theta(std::atof(option_value.c_str()));
      break;
    case GAMMA_F:
      target().gamma_f(std::atof(option_value.c_str()));
      break;
    case F_MIN: {
      if( option_value == "UNBOUNDED" )
        target().f_min(target_t::F_MIN_UNBOUNDED);
      else
        target().f_min(std::atof(option_value.c_str()));
      break;
    }
    case GAMMA_ALPHA:
      target().gamma_alpha(std::atof(option_value.c_str()));
      break;
    case DELTA:
      target().delta(std::atof(option_value.c_str()));
      break;
    case S_THETA:
      target().s_theta(std::atof(option_value.c_str()));
      break;
    case S_F:
      target().s_f(std::atof(option_value.c_str()));
      break;
    case THETA_SMALL_FACT:
      target().theta_small_fact(std::atof(option_value.c_str()));
      break;
    case THETA_MAX:
      target().theta_max(std::atof(option_value.c_str()));
      break;
    case ETA_F:
      target().eta_f(std::atof(option_value.c_str()));
      break;
    case BACK_TRACK_FRAC:
      target().back_track_frac(std::atof(option_value.c_str()));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
    }
  }

}	// end namespace MoochoPack 
