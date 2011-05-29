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

#ifndef DECOMPOSITION_SYSTEM_TESTER_SET_OPTIONS_H
#define DECOMPOSITION_SYSTEM_TESTER_SET_OPTIONS_H

#include "ConstrainedOptPack_DecompositionSystemTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for DecompositionSystemTester from an
  * OptionsFromStream object.
  *
  * The default options group name is DecompositionSystemTester.
  *
  * The options group is:
  *
  \verbatim

    options_group DecompositionSystemTester {
        print_tests = PRINT_NONE;
    *    print_tests = PRINT_BASIC;
    *    print_tests = PRINT_MORE;
    *    print_tests = PRINT_ALL;
    *    dump_all = true;
        dump_all = false;
        throw_exception = true;
    *    throw_exception = false;
        num_random_tests = 1;
        mult_warning_tol = 1e-14;
        mult_error_tol   = 1e-10;
        solve_warning_tol = 1e-14;
        solve_error_tol   = 1e-10;
    }
  \endverbatim
  */
class DecompositionSystemTesterSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      DecompositionSystemTester >
{
public:

  /** \brief . */
  DecompositionSystemTesterSetOptions(
      DecompositionSystemTester* target = 0
    , const char opt_grp_name[] = "DecompositionSystemTester" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class DecompositionSystemTesterSetOptions

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_TESTER_SET_OPTIONS_H
