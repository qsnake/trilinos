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

#include "Moocho_ConfigDefs.hpp"

#ifdef HAVE_MOOCHO_MA28

#include "AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 8;

  enum local_EOptions {
    ESTIMATED_FILLIN_RATIO
      ,U
    ,GROW
    ,TOL
    ,NSRCH
    ,LBIG
    ,PRINT_MA28_OUTPUTS
    ,OUTPUT_FILE_NAME
  };

  const char* local_SOptions[local_num_options]	= {
    "estimated_fillin_ratio"
      ,"u"
    ,"grow"
    ,"tol"
    ,"nsrch"
    ,"lbig"
    ,"print_ma28_outputs"
    ,"output_file_name"
  };

}

namespace AbstractLinAlgPack {

DirectSparseSolverMA28SetOptions::DirectSparseSolverMA28SetOptions(
      DirectSparseSolverMA28* qp_solver )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        "DirectSparseSolverMA28", local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      DirectSparseSolverMA28 >( qp_solver )
{}

void DirectSparseSolverMA28SetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  switch( (local_EOptions)option_num ) {
    case ESTIMATED_FILLIN_RATIO: {
      target().estimated_fillin_ratio( ::atof( option_value.c_str() ) );
      break;
    }
    case U: {
      target().u( ::atof( option_value.c_str() ) );
      break;
    }
    case GROW: {
      target().grow( StringToBool( "grow", option_value.c_str() ) );
      break;
    }
    case TOL: {
      target().tol( ::atof( option_value.c_str() ) );
      break;
    }
    case NSRCH: {
      target().nsrch( ::atoi( option_value.c_str() ) );
      break;
    }
    case LBIG: {
      target().lbig( StringToBool( "lbig", option_value.c_str() ) );
      break;
    }
    case PRINT_MA28_OUTPUTS: {
      target().print_ma28_outputs( StringToBool( "grow", option_value.c_str() ) );
      break;
    }
    case OUTPUT_FILE_NAME: {
      if( option_value == "NONE" )
        target().output_file_name( "" );
      else
        target().output_file_name( option_value );
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace AbstractLinAlgPack 

#endif // HAVE_MOOCHO_MA28
