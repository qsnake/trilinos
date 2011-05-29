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

#ifndef EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H
#define EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H

#include "MoochoPack_EvalNewPointStd_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for EvalNewPointStd_Step from an \c OptionsFromStream object.
 *
 * The default options group name is EvalNewPointStd.
 *
 * The options group is:
 *
 \verbatim
    options_group EvalNewPointStd {
        fd_deriv_testing = FD_DEFAULT;
        decomp_sys_teting = DST_DEFAULT;
        decomp_sys_teting_print_level = DSPL_USE_GLOBAL;
    }
 \verbatim
 *
 * <ul>
 * <li> <b>fd_deriv_testing</b>: Determines if finite differerece testing of the 
 *      derivatives of the Gc and Gf.  See the class \c EvalNewPointStd_Step
 *      and its printed algorithm for more details.
 *      <ul>
 *      <li> <b>FD_DEFAULT</b>: The global flag check_results determines
 *           if the tests are performed.
 *      <li> <b>FD_TEST</b>: The tests are performed reguardless the
 *           value of check_results
 *      <li> <b>FD_NO_TEST</b>: The tests are not performed reguardless the
 *           value of check_results
 *      </ul>
 * <li> <b>decomp_sys_testing</b>: Determines if the range/null decomposition of
 *      Gc and Gh is performed.  See the class \c EvalNewPointStd_Step
 *      and its printed algorithm for more details.
 *      <ul>
 *      <li> <b>DST_DEFAULT</b>: The global flag check_results determines
 *           if the tests are performed.
 *      <li> <b>DST_TEST</b>: The tests are performed reguardless the
 *           value of check_results
 *      <li> <b>DST_NO_TEST</b>: The tests are not performed reguardless the
 *           value of check_results
 *      </ul>
 * </ul>
 */
class EvalNewPointStd_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      EvalNewPointStd_Step >
{
public:

  /** \brief . */
  EvalNewPointStd_StepSetOptions(
      EvalNewPointStd_Step* target = 0
    , const char opt_grp_name[] = "EvalNewPointStd" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class EvalNewPointStd_StepSetOptions

}	// end namespace MoochoPack

#endif	// EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H
