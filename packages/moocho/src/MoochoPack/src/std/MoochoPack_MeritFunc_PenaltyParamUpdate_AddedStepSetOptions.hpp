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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for MeritFunc_PenaltyParamUpdate_AddedStep from a
 * OptionsFromStream object.
 *
 * The options group is:
 *
 \verbatim
  options_group MeritFuncPenaltyParamUpdate {
      small_mu     = 1e-6;
      min_mu_ratio = 1e-8;
      mult_factor  = 1e-4;
      kkt_near_sol = 1.0;
  }
 \endverbatim
 *
 * <ul>
 *	<li>[small_mu] The smallest mu allows when away from the soltion.<br>
 *		Example: small_mu = 1e-6;
 *	<li>[min_mu_ratio] This bounds the smallest mu(i) as:<br>
 *		min(mu(i))/max(mu(i)) >= min_mu_ratio.<br>
 *		Example: min_mu_ratio = 1e-4;
 *	<li>[mult_factor] Multiplicative factor for mu(j) = (1.0+mult_factor) * abs(lambda(j)).<br>
 *		Example: mult_factor = 1e-4;
 *	<li>[kkt_near_sol] When the total kkt_error is below kkt_near_sol a safer
 *		penalty update will be used.<br>
 *		Example: kkt_near_sol = 1.0;
 * </ul>
 */
class MeritFunc_PenaltyParamUpdate_AddedStepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      MeritFunc_PenaltyParamUpdate_AddedStep >
{
public:

  /** \brief . */
  MeritFunc_PenaltyParamUpdate_AddedStepSetOptions(
    MeritFunc_PenaltyParamUpdate_AddedStep* target = 0 );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class MeritFunc_PenaltyParamUpdate_AddedStepSetOptions

}	// end namespace MoochoPack

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H
