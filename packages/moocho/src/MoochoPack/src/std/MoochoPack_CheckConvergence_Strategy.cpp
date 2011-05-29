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

#include "MoochoPack_CheckConvergence_Strategy.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

namespace MoochoPack {

///*******************************************
//  CheckConvergence_Strategy
///*******************************************


CheckConvergence_Strategy::CheckConvergence_Strategy(
  EOptErrorCheck opt_error_check,
  EScaleKKTErrorBy scale_opt_error_by,
  EScaleKKTErrorBy scale_feas_error_by,
  EScaleKKTErrorBy scale_comp_error_by,
  bool scale_opt_error_by_Gf
  )
  :
  opt_error_check_(opt_error_check),
  scale_opt_error_by_(scale_opt_error_by),
  scale_feas_error_by_(scale_feas_error_by),
  scale_comp_error_by_(scale_comp_error_by),
  scale_opt_error_by_Gf_(scale_opt_error_by_Gf) 
  {}


///*******************************************
//  CheckConvergence_StrategySetOptions
///*******************************************

// Define the options
namespace {

const int local_num_options = 4;

enum local_EOptions 
  {
  SCALE_OPT_ERROR_BY,
  SCALE_FEAS_ERROR_BY,
  SCALE_COMP_ERROR_BY,
  SCALE_OPT_ERROR_BY_GF
  };

const char* local_SOptions[local_num_options] = 
  {
  "scale_opt_error_by",
  "scale_feas_error_by",
  "scale_comp_error_by",
  "scale_opt_error_by_Gf",
  };

} // end namespace

CheckConvergence_StrategySetOptions::CheckConvergence_StrategySetOptions(
  CheckConvergence_Strategy* target,
  const char opt_grp_name[] )
  :  
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase<
    CheckConvergence_Strategy >( target )
  {}


void CheckConvergence_StrategySetOptions::setOption(
  int option_num,
  const std::string& option_value )
  {
  using OptionsFromStreamPack::StringToBool;

  typedef CheckConvergence_Strategy target_t;
  switch( (local_EOptions)option_num ) 
    {
    case SCALE_OPT_ERROR_BY:
    case SCALE_FEAS_ERROR_BY:
    case SCALE_COMP_ERROR_BY:
      {
      const std::string &option = option_value.c_str();
      CheckConvergence_Strategy::EScaleKKTErrorBy scale_by = target_t::SCALE_BY_ONE;

      if( option == "SCALE_BY_ONE" )
        { scale_by = target_t::SCALE_BY_ONE; }
      else if( option == "SCALE_BY_NORM_2_X" )
        { scale_by = target_t::SCALE_BY_NORM_2_X; }
      else if( option == "SCALE_BY_NORM_INF_X" )
        { scale_by = target_t::SCALE_BY_NORM_INF_X; }
      else
        {
        throw std::invalid_argument( "Error, incorrect value for "
                       "\"scale_kkt_error_by\".  Only the options "
                       "SCALE_BY_ONE, SCALE_BY_NORM_2_X, and SCALE_BY_NORM_INF_X "
                       "are available" );
        }


      if ((local_EOptions) option_num == SCALE_OPT_ERROR_BY)
        {
        target().scale_opt_error_by(scale_by);
        }
      else if ((local_EOptions) option_num == SCALE_FEAS_ERROR_BY)
        {
        target().scale_feas_error_by(scale_by);
        }
      else if ((local_EOptions) option_num == SCALE_COMP_ERROR_BY)
        {
        target().scale_comp_error_by(scale_by);
        }
      else
        {
        TEST_FOR_EXCEPTION( true,
                 std::logic_error,
                 "Unaccounted for option_num in CheckConvergence_Strategy.cpp"
          );
        }

      break;
      }
    case SCALE_OPT_ERROR_BY_GF: 
      {
      target().scale_opt_error_by_Gf(
        StringToBool( "scale_opt_error_by_Gf", option_value.c_str() ) );
      break;
      }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
    }
  }

} // end namespace MoochoPack
