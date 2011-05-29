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

#include "ConstrainedOptPack_VariableBoundsTesterSetOptions.hpp"

// Define the options
namespace {

  const int local_num_options = 2;

  enum local_EOptions {
    WARNING_TOL
    ,ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
      "warning_tol"
      ,"error_tol"
  };

}

namespace ConstrainedOptPack {

VariableBoundsTesterSetOptions::VariableBoundsTesterSetOptions(
        VariableBoundsTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      VariableBoundsTester >( target )
{}

void VariableBoundsTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef VariableBoundsTester target_t;
  switch( (local_EOptions)option_num ) {
      case WARNING_TOL:
      target().warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case ERROR_TOL:
      target().error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack
