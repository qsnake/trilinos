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

#include "NLPInterfacePack_NLPTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 2;

  enum local_EOptions {
    PRINT_ALL
    ,TEST_FOR_EXCEPTION
  };

  const char* local_SOptions[local_num_options]	= {
    "print_all"
    ,"throw_exception"
  };

}

namespace NLPInterfacePack {

NLPTesterSetOptions::NLPTesterSetOptions(
        NLPTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      NLPTester >( target )
{}

void NLPTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  namespace ofsp = OptionsFromStreamPack;
  using ofsp::StringToBool;

  typedef NLPTester target_t;
  switch( (local_EOptions)option_num ) {
    case PRINT_ALL:
      target().print_all(
        StringToBool( "print_all", option_value.c_str() )
        );
      break;
    case TEST_FOR_EXCEPTION:
      target().throw_exception(
        StringToBool( "throw_exception", option_value.c_str() )
        );
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace NLPInterfacePack
