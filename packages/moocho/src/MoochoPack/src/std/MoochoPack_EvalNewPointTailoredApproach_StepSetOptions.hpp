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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H

#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for EvalNewPointTailoredApproach_Step from an
  * OptionsFromStream object.
  *
  * The default options group name is EvalNewPointTailoredApproach.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group EvalNewPointTailoredApproach {
    fd_deriv_testing   = FD_DEFAULT;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[fd_deriv_testing] Determines if finite differerece testing of the 
  *		derivatives of the Gc and Gf.  See the class \Ref{EvalNewPointTailoredApproach_Step}
  *		and its printed algorithm for more details).
  *		\begin{description}
  *		\item[FD_DEFAULT]			The global flag check_results determines
  *									if the tests are performed.
  *		\item[FD_TEST]				The tests are performed reguardless the
  *									value of check_results
  *		\item[FD_NO_TEST]			The tests are not performed reguardless the
  *									value of check_results
  *		\end{description}
  */
class EvalNewPointTailoredApproach_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      EvalNewPointTailoredApproach_Step >
{
public:

  /** \brief . */
  EvalNewPointTailoredApproach_StepSetOptions(
      EvalNewPointTailoredApproach_Step* target = 0
    , const char opt_grp_name[] = "EvalNewPointTailoredApproach" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class EvalNewPointTailoredApproach_StepSetOptions

}	// end namespace MoochoPack

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_STD_SET_OPTIONS_H
