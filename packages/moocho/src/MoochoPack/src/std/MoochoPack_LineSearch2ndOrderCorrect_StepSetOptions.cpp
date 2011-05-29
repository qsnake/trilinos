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

#include "MoochoPack_LineSearch2ndOrderCorrect_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 8;

  enum local_EOptions {
      NEWTON_OLEVEL,
    CONSTR_NORM_THRESHOLD,
    CONSTR_INCR_RATIO,
    AFTER_K_ITER,
    FORCED_CONSTR_REDUCTION,
    FORCED_REDUCT_RATIO,
    MAX_STEP_RATIO,
    MAX_NEWTON_ITER
  };

  const char* local_SOptions[local_num_options]	= {
      "newton_olevel",
    "constr_norm_threshold",
    "constr_incr_ratio",
    "after_k_iter",
    "forced_constr_reduction",
    "forced_reduct_ratio",
    "max_step_ratio",
    "max_newton_iter"
  };

}

namespace MoochoPack {

LineSearch2ndOrderCorrect_StepSetOptions::LineSearch2ndOrderCorrect_StepSetOptions(
        LineSearch2ndOrderCorrect_Step* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      LineSearch2ndOrderCorrect_Step >( target )
{}

void LineSearch2ndOrderCorrect_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef LineSearch2ndOrderCorrect_Step target_t;
  switch( (local_EOptions)option_num ) {
      case NEWTON_OLEVEL:
    {
      const std::string &option = option_value.c_str();
      if( option == "PRINT_USE_DEFAULT" )
        target().newton_olevel( target_t::PRINT_USE_DEFAULT );
      else if( option == "PRINT_NOTHING" )
        target().newton_olevel( target_t::PRINT_NEWTON_NOTHING );
      else if( option == "PRINT_SUMMARY_INFO" )
        target().newton_olevel( target_t::PRINT_NEWTON_SUMMARY_INFO );
      else if( option == "PRINT_STEPS" )
        target().newton_olevel( target_t::PRINT_NEWTON_STEPS );
      else if( option == "PRINT_VECTORS" )
        target().newton_olevel( target_t::PRINT_NEWTON_VECTORS );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"newton_olevel\"." );
      break;
    }
      case CONSTR_NORM_THRESHOLD:
      target().constr_norm_threshold(::fabs(::atof(option_value.c_str())));
      break;
      case CONSTR_INCR_RATIO:
      target().constr_incr_ratio(::fabs(::atof(option_value.c_str())));
      break;
    case AFTER_K_ITER:
      target().after_k_iter(::abs(::atoi(option_value.c_str())));
      break;
    case FORCED_CONSTR_REDUCTION:
    {
      const std::string &option = option_value.c_str();
      if( option == "LESS_X_D" )
        target().forced_constr_reduction(target_t::CONSTR_LESS_X_D );
      else if( option == "LESS_X" )
        target().forced_constr_reduction( target_t::CONSTR_LESS_X );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"forced_constr_reduction\"." );
      break;
    }
      case FORCED_REDUCT_RATIO:
      target().forced_reduct_ratio(::fabs(::atof(option_value.c_str())));
      break;
    case MAX_STEP_RATIO:
      target().max_step_ratio(::fabs(::atof(option_value.c_str())));
      break;
    case MAX_NEWTON_ITER:
      target().max_newton_iter(::abs(::atoi(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 

#endif // 0
