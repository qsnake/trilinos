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

#include "IterationPack_AlgorithmSetOptions.hpp"
#include "Teuchos_TestForException.hpp"

// Define the options
namespace {

const int local_num_options = 1;

enum local_EOptions {
  INTERRUPT_FILE_NAME
};

const char* local_SOptions[local_num_options]	= {
  "interrupt_file_name"
};

inline
std::string remove_quotes( const std::string &option_name, const std::string &str )
{
  TEST_FOR_EXCEPTION(
    str[0]!='\"' || str[str.length()-1]!='\"', std::logic_error
    ,"Error, the option \'" << option_name << "\' must have a single set of quotes around it!"
    );
  if(str.length()==2)
    return std::string("");
  return str.substr(1,str.length()-2);
}

} // namespace

namespace IterationPack {

AlgorithmSetOptions::AlgorithmSetOptions(
        Algorithm* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      Algorithm >( target )
{}

void AlgorithmSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef Algorithm target_t;
  switch( (local_EOptions)option_num ) {
    case INTERRUPT_FILE_NAME :
      target().interrupt_file_name(remove_quotes("interrupt_file_name",option_value));
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace IterationPack 
