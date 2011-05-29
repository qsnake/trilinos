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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H
#define REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_ReducedHessianSecantUpdateLPBFGS_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for ReducedHessianSecantUpdateBFGSProjected_Strategy
  * from a OptionsFromStream object.
  *
  * The options group is (with the default name):
  *
  \begin{verbatim}
  options_group ReducedHessianSecantUpdateLPBFGS {
    min_num_updates_proj_start   = 0;      *** (+int)
    max_num_updates_proj_start   = 999999; *** (+int)
    num_superbasics_switch_dense = 500;    *** (+int)
    num_add_recent_updates       = 10;     *** (+int)
    }
  \end{verbatim}
  *
  * \begin{description}
  *	\item ToDo : Finish
  *	\end{description}
  */
class ReducedHessianSecantUpdateLPBFGS_StrategySetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
  , public OptionsFromStreamPack::SetOptionsToTargetBase<
        ReducedHessianSecantUpdateLPBFGS_Strategy >
  {
public:

  /** \brief . */
  ReducedHessianSecantUpdateLPBFGS_StrategySetOptions(
    ReducedHessianSecantUpdateLPBFGS_Strategy* target = 0
    , const char opt_grp_name[] = "ReducedHessianSecantUpdateLPBFGS" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateLPBFGS_StrategySetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_SET_OPTIONS_H
