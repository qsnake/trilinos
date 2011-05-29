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

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_FeasibilityStepReducedStd_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for FeasibilityStepReducedStd_Strategy from an
  * OptionsFromStream object.
  *
  * The default options group name is IndepDirecWithBoundsStd.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group FeasibilityStepReducedStd_Strategy {
  *    qp_objective = OBJ_MIN_FULL_STEP;
  *    qp_objective = OBJ_MIN_NULL_SPACE_STEP;
  *    qp_objective = OBJ_RSQP;
  *    qp_testing   = QP_TEST_DEFAULT;
  *    qp_testing   = QP_TEST;
  *    qp_testing   = QP_NO_TEST;
  }
  \end{verbatim}
  */
class FeasibilityStepReducedStd_StrategySetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      FeasibilityStepReducedStd_Strategy >
{
public:

  /** \brief . */
  FeasibilityStepReducedStd_StrategySetOptions(
      FeasibilityStepReducedStd_Strategy* target = 0
    , const char opt_grp_name[] = "FeasibilityStepReducedStd" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class FeasibilityStepReducedStd_StrategySetOptions

}	// end namespace MoochoPack

#endif	// FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
