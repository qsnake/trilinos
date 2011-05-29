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

#ifndef VECTOR_SPACE_TESTER_SET_OPTIONS_H
#define VECTOR_SPACE_TESTER_SET_OPTIONS_H

#include "AbstractLinAlgPack_VectorSpaceTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace AbstractLinAlgPack {

/** \brief Set options for VectorSpaceTester from an
  * OptionsFromStream object.
  *
  * The default options group name is VectorSpaceTester.
  *
  * The options group is:
  *
  \verbatim

    options_group VectorSpaceTester {
    *    print_all_tests = true;
        print_all_tests = false;
    *    print_vectors = true;
        print_vectors = false;
        throw_exception = true;
    *    throw_exception = false;
        num_random_tests = 4;
        warning_tol = 1e-14;
        error_tol   = 1e-10;
    }
  \endverbatim
  */
class VectorSpaceTesterSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      VectorSpaceTester >
{
public:

  /** \brief . */
  VectorSpaceTesterSetOptions(
      VectorSpaceTester* target = 0
    , const char opt_grp_name[] = "VectorSpaceTester" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class VectorSpaceTesterSetOptions

}	// end namespace AbstractLinAlgPack

#endif	// VECTOR_SPACE_TESTER_SET_OPTIONS_H
