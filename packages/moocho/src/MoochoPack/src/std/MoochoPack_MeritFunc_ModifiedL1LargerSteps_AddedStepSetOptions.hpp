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

#ifndef MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H
#define MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H

#include "MoochoPack_MeritFunc_ModifiedL1LargerSteps_AddedStep.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for MeritFunc_ModifiedL1LargerSteps_AddedStep from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group MeritFuncModifiedL1LargerSteps {
    after_k_iter                = 3;
    obj_increase_threshold      = 1e-3;
    max_pos_penalty_increase    = 1.0;
    pos_to_neg_penalty_increase = 1.0;
    incr_mult_factor            = 1e-4;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[after_k_iter] ToDo : Finish.
  *		Example: after_k_iter = 4;
  *	\item[obj_increase_threshold] ToDo : Finish.
  *		Example: obj_increase_threshold = 1e-4;
  *	\item[max_pos_penalty_increase] ToDo : Finish.
  *		Example: max_pos_penalty_increase = 1.0;
  *	\item[pos_to_neg_penalty_increase] ToDo : Finish.
  *		Example: pos_to_neg_penalty_increase = 1.0;
  *	\item[incr_mult_factor] ToDo : Finish.
  *		Example: incr_mult_factor = 1e-4;
  *	\end{description}
  */
class MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      MeritFunc_ModifiedL1LargerSteps_AddedStep >
{
public:

  /** \brief . */
  MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions(
    MeritFunc_ModifiedL1LargerSteps_AddedStep* target = 0 );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions

}	// end namespace MoochoPack

#endif	// MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H
