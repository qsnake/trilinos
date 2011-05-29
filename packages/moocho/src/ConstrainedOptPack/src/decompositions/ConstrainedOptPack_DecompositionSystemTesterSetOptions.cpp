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

#include "ConstrainedOptPack_DecompositionSystemTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

  const int local_num_options = 8;

  enum local_EOptions {
    PRINT_TESTS
    ,DUMP_ALL
    ,TEST_FOR_EXCEPTION
    ,NUM_RANDOM_TESTS
      ,MULT_WARNING_TOL
      ,MULT_ERROR_TOL
      ,SOLVE_WARNING_TOL
      ,SOLVE_ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
    "print_tests"
    ,"dump_all"
    ,"throw_exception"
    ,"num_random_tests"
      ,"mult_warning_tol"
      ,"mult_error_tol"
      ,"solve_warning_tol"
      ,"solve_error_tol"
  };

}

namespace ConstrainedOptPack {

DecompositionSystemTesterSetOptions::DecompositionSystemTesterSetOptions(
        DecompositionSystemTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      DecompositionSystemTester >( target )
{}

void DecompositionSystemTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  typedef DecompositionSystemTester target_t;
  switch( (local_EOptions)option_num ) {
    case PRINT_TESTS:
    {
      const std::string &option = option_value.c_str();
      if( option == "PRINT_NONE" )
        target().print_tests( target_t::PRINT_NONE );
        else if( option == "PRINT_BASIC" )
        target().print_tests( target_t::PRINT_BASIC );
        else if( option == "PRINT_MORE" )
        target().print_tests( target_t::PRINT_MORE );
        else if( option == "PRINT_ALL" )
        target().print_tests( target_t::PRINT_ALL );
      else
        TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"Error, incorrect value for "
          "\"print_tests\".  Only the options "
          "PRINT_NONE, PRINT_BASIS, PRINT_MORE and PRINT_ALL are allowed" );
      break;
    }
    case DUMP_ALL:
      target().dump_all(
        StringToBool( "dump_all", option_value.c_str() )
        );
      break;
    case TEST_FOR_EXCEPTION:
      target().throw_exception(
        StringToBool( "throw_exception", option_value.c_str() )
        );
      break;
      case NUM_RANDOM_TESTS:
      target().num_random_tests(std::abs(std::atoi(option_value.c_str())));
      break;
      case MULT_WARNING_TOL:
      target().mult_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case MULT_ERROR_TOL:
      target().mult_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case SOLVE_WARNING_TOL:
      target().solve_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case SOLVE_ERROR_TOL:
      target().solve_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack
