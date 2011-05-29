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

#ifndef LINE_SEARCH_WATCH_DOG_STEP_H
#define LINE_SEARCH_WATCH_DOG_STEP_H

#include "rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "MiStandardAggregationMacros.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implements watchdog line search.
  *
  * The watchdog proceedure will only be considered when :
  * opt_kkt_err_k <= opt_kkt_err_threshold && feas_kkt_err_k <= feas_kkt_err_threshold
  * (see step listing).  The default behavior is to never use the watchdog procedure.
  */
class LineSearchWatchDog_Step : public LineSearch_Step {
public:

  /// <<std comp>> members for direct_line_search
  STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search);

  /// <<std comp>> members for merit_func
  STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func);

  /** \brief <<std member comp>> members for the armijo fractional reduction parameter.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta );

  /** \brief <<std member comp>> members for the threshold for opt_kkt_err before
    * the watchdog procedure should kick-in.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_kkt_err_threshold );
  
  /** \brief <<std member comp>> members for the threshold for feas_kkt_err before
    * the watchdog procedure should kick-in.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_kkt_err_threshold );

  /** \brief . */
  LineSearchWatchDog_Step(
        const direct_line_search_ptr_t&	direct_line_search		= 0
      , const merit_func_ptr_t&			merit_func				= 0
      , value_type						eta						= 1e-4
      , value_type						opt_kkt_err_threshold 	= 1e-1
      , value_type						feas_kkt_err_threshold 	= 1e-3
      );

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
  int					watch_k_;	// >= 0 means that we are using the watchdog.
  DVector				xo_;
  value_type			fo_;
  value_type			nrm_co_;
  DVector				do_;
  value_type			phio_;
  value_type			Dphio_;
  value_type			phiop1_;
  
};	// end class LineSearchWatchDog_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_WATCH_DOG_STEP_H
