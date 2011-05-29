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

#include "ConstrainedOptPack_DirectLineSearchArmQuad_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 5;

  enum local_EOptions {
    SLOPE_FRAC,
    MIN_FRAC_STEP,
    MAX_FRAC_STEP,
    MAX_LS_ITER,
    MAX_OUT_LS_ITER
  };

  const char* local_SOptions[local_num_options]	= {
    "slope_frac",
    "min_frac_step",
    "max_frac_step",
    "max_ls_iter",
    "max_out_ls_iter"
  };

}

namespace ConstrainedOptPack {

DirectLineSearchArmQuad_StrategySetOptions::DirectLineSearchArmQuad_StrategySetOptions(
        DirectLineSearchArmQuad_Strategy* qp_solver
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      DirectLineSearchArmQuad_Strategy >( qp_solver )
{}

void DirectLineSearchArmQuad_StrategySetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  switch( (local_EOptions)option_num ) {
    case SLOPE_FRAC: {
      target().eta( std::atof( option_value.c_str() ) );
      break;
    }
    case MIN_FRAC_STEP: {
      target().min_frac( std::atof( option_value.c_str() ) );
      break;
    }
    case MAX_FRAC_STEP: {
      target().max_frac( std::atof( option_value.c_str() ) );
      break;
    }
    case MAX_LS_ITER: {
      target().set_max_iter( std::atof( option_value.c_str() ) );
      break;
    }
    case MAX_OUT_LS_ITER: {
      target().max_out_iter( StringToBool( "max_out_ls_iter", option_value.c_str() ) );
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack 
