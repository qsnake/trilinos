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

#ifndef RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H
#define RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H

#include "MoochoPack_NLPSolverClientInterface.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for NLPSolverClientInterface from an \c OptionsFromStream object.
 *
 * The default options group name is NLPSolverClientInterface.
 *
 * The options group is:
 *
 \verbatim
  options_group NLPSolverClientInterface {
        max_iter = ?;
        max_run_time = ?;  *** In minutes
        opt_tol = ?;
        feas_tol = ?;
        step_tol = ?;
    journal_output_level = ?;
    journal_print_digits = ?;
    check_results = ?;
    calc_conditioning = ?
  }
 \endverbatim
 *
 * See the class \c NLPSolverClientInterface for a description of these
 * parameters.
 */
class NLPSolverClientInterfaceSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      NLPSolverClientInterface >
{
public:

  /** \brief . */
  NLPSolverClientInterfaceSetOptions(
      NLPSolverClientInterface* target = 0
    , const char opt_grp_name[] = "NLPSolverClientInterface" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class NLPSolverClientInterfaceSetOptions

}	// end namespace MoochoPack

#endif	// RSQP_SOLVER_CLIENT_INTERFACE_SET_OPTIONS_H
