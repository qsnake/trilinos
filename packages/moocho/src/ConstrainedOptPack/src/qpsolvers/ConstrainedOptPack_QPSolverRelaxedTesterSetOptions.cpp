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

#include "ConstrainedOptPack_QPSolverRelaxedTesterSetOptions.hpp"

// Define the options
namespace {

  const int local_num_options = 6;

  enum local_EOptions {
    OPT_WARNING_TOL
    ,OPT_ERROR_TOL
    ,FEAS_WARNING_TOL
    ,FEAS_ERROR_TOL
    ,COMP_WARNING_TOL
    ,COMP_ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
      "opt_warning_tol"
      ,"opt_error_tol"
      ,"feas_warning_tol"
      ,"feas_error_tol"
      ,"comp_warning_tol"
      ,"comp_error_tol"
  };

}

namespace ConstrainedOptPack {

QPSolverRelaxedTesterSetOptions::QPSolverRelaxedTesterSetOptions(
        QPSolverRelaxedTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      QPSolverRelaxedTester >( target )
{}

void QPSolverRelaxedTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef QPSolverRelaxedTester target_t;
  switch( (local_EOptions)option_num ) {
      case OPT_WARNING_TOL:
      target().opt_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case OPT_ERROR_TOL:
      target().opt_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case FEAS_WARNING_TOL:
      target().feas_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case FEAS_ERROR_TOL:
      target().feas_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case COMP_WARNING_TOL:
      target().comp_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case COMP_ERROR_TOL:
      target().comp_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack
