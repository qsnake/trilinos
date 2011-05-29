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

#ifndef DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_SET_OPTIONS_H
#define DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_SET_OPTIONS_H

#include "ConstrainedOptPack_DirectLineSearchArmQuad_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for DirectLineSearchArmQuad_Strategy from a
  * OptionsFromStream object.
  *
  * The default options group name is DirectLineSearchArmQuad.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group DirectLineSearchArmQuad {
    slope_frac = ?;
    min_frac_step = ?:
    max_frac_step = ?;
    max_ls_iter = ?;
    max_out_ls_iter = ?;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[slope_frac] Fraction of intial decent slope of merit
  *		function required by armijo test.  Mapps to eta in linesearch
  *		algorithm.  Must be in the range [0, 1].
  *		Example: slope_frac = 1.0e-4;
  *	\item[min_frac_step] Minimum fractional change in the step length
  *		per linesearch iteration.  Mapps to min_frac in the linesearch
  *		algorithm.  Must be in the range (0, 1].
  *		Example: min_frac_step = 0.1;
  *	\item[max_frac_step] Maximum fractional change in the step length
  *		per linesearch iteration.  Mapps to max_frac in the linesearch
  *		algorithm.  Must be in the range (0, 1].  Note that
  *		max_frac_step must be greater than min_frac_step.
  *		Example: min_frac_step = 0.5;
  *	\item[max_ls_iter] The maximum number of linesearch iterations
  *		to take before giving up and declaring a line search failure.
  *		Mapps to max_ls_iter in linesearch algorithm.
  *		Example: max_ls_iter = 20.
  *	\item[max_out_ls_iter] A flag to max out on line search iterations.
  *   Mostly just used for debugging, not very useful in general.
  *	\end{description}
  */
class DirectLineSearchArmQuad_StrategySetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      DirectLineSearchArmQuad_Strategy >
{
public:

  /** \brief . */
  DirectLineSearchArmQuad_StrategySetOptions(
      DirectLineSearchArmQuad_Strategy* qp_solver = 0
    , const char opt_grp_name[] = "DirectLineSearchArmQuad" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class DirectLineSearchArmQuad_StrategySetOptions

}	// end namespace ConstrainedOptPack

#endif	// DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_SET_OPTIONS_H
