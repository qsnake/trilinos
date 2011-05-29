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

#ifndef LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
#define LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H

#include "MoochoPack_LineSearchWatchDog_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for LineSearchWatchDog_Step from a OptionsFromStream
  * object.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group LineSearchWatchDog {
    opt_kkt_err_threshold	= 1e-3; *** (+dbl)
    feas_kkt_err_threshold	= 1e-3; *** (+dbl)
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[opt_kkt_err_threshold] ToDo : Finish.
  *		Example: opt_kkt_err_threshold = 1e-1;
  *	\item[feas_kkt_err_threshold] ToDo : Finish.
  *		Example: feas_kkt_err_threshold = 1e-2;
  *	\end{description}
  */
class LineSearchWatchDog_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      LineSearchWatchDog_Step >
{
public:

  /** \brief . */
  LineSearchWatchDog_StepSetOptions(
    LineSearchWatchDog_Step* target = 0 );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class LineSearchWatchDog_StepSetOptions

}	// end namespace MoochoPack

#endif	// LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
