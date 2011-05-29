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

#include "NLPInterfacePack_NLPDirectTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 8;

  enum local_EOptions {
    GF_TESTING_METHOD
    ,GF_WARNING_TOL
    ,GF_ERROR_TOL
    ,GC_TESTING_METHOD
    ,GC_WARNING_TOL
    ,GC_ERROR_TOL
    ,NUM_FD_DIRECTIONS
    ,DUMP_ALL
  };

  const char* local_SOptions[local_num_options]	= {
      "Gf_testing_method"
      ,"Gf_warning_tol"
      ,"Gf_error_tol"
      ,"Gc_testing_method"
      ,"Gc_warning_tol"
      ,"Gc_error_tol"
      ,"num_fd_directions"
      ,"dump_all"
  };

}

namespace NLPInterfacePack {

NLPDirectTesterSetOptions::NLPDirectTesterSetOptions(
  NLPDirectTester* target
  ,const char opt_grp_name[]
  )
  :OptionsFromStreamPack::SetOptionsFromStreamNode(opt_grp_name,local_num_options,local_SOptions)
  ,OptionsFromStreamPack::SetOptionsToTargetBase<NLPDirectTester>( target )
{}

void NLPDirectTesterSetOptions::setOption(
  int option_num, const std::string& option_value
  )
{
  namespace ofsp = OptionsFromStreamPack;
  using ofsp::StringToBool;
  typedef NLPDirectTester target_t;
  switch( (local_EOptions)option_num ) {
    case GF_TESTING_METHOD:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_COMPUTE_ALL" )
        target().Gf_testing_method( target_t::FD_COMPUTE_ALL );
      else if( option == "FD_DIRECTIONAL" )
        target().Gf_testing_method( target_t::FD_DIRECTIONAL );
      else
        throw std::invalid_argument( "Error, incorrect value for "
                                     "\"Gf_testing_method\".  Only the options "
                                     "FD_COMPUTE_ALL and FD_DIRECTIONAL are available" );
      break;
    }
    case GF_WARNING_TOL:
      target().Gf_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case GF_ERROR_TOL:
      target().Gf_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case GC_TESTING_METHOD:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_COMPUTE_ALL" )
        target().Gc_testing_method( target_t::FD_COMPUTE_ALL );
      else if( option == "FD_DIRECTIONAL" )
        target().Gc_testing_method( target_t::FD_DIRECTIONAL );
      else
        throw std::invalid_argument( "Error, incorrect value for "
                                     "\"Gc_testing_method\".  Only the options "
                                     "FD_COMPUTE_ALL and FD_DIRECTIONAL are available" );
      break;
    }
    case GC_WARNING_TOL:
      target().Gc_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case GC_ERROR_TOL:
      target().Gc_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case NUM_FD_DIRECTIONS:
      target().num_fd_directions(std::abs(std::atoi(option_value.c_str())));
      break;
    case DUMP_ALL:
      target().dump_all(StringToBool("dump_all",option_value.c_str()));
      break;
    default:
      TEST_FOR_EXCEPT(true); // Local error only?
  }
}

}	// end namespace NLPInterfacePack
