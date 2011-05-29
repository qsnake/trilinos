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

#include "AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 6;

  enum local_EOptions {
    PRINT_ALL_TESTS
    ,PRINT_VECTORS
    ,TEST_FOR_EXCEPTION
    ,NUM_RANDOM_TESTS
      ,WARNING_TOL
      ,ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
    "print_all_tests"
    ,"print_vectors"
    ,"throw_exception"
    ,"num_random_tests"
      ,"warning_tol"
      ,"error_tol"
  };

}

namespace AbstractLinAlgPack {

VectorSpaceTesterSetOptions::VectorSpaceTesterSetOptions(
        VectorSpaceTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      VectorSpaceTester >( target )
{}

void VectorSpaceTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  typedef VectorSpaceTester target_t;
  switch( (local_EOptions)option_num ) {
    case PRINT_ALL_TESTS:
      target().print_all_tests(
        StringToBool( "print_all_tests", option_value.c_str() )
        );
      break;
    case PRINT_VECTORS:
      target().print_vectors(
        StringToBool( "print_vectors", option_value.c_str() )
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
      case WARNING_TOL:
      target().warning_tol(std::abs(std::atof(option_value.c_str())));
      break;
      case ERROR_TOL:
      target().error_tol(std::abs(std::atof(option_value.c_str())));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace AbstractLinAlgPack
