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

#ifndef LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
#define LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H

#include <limits>

#include "rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "MiStandardAggregationMacros.h"

namespace MoochoPack {

/** \brief Changes from a line search step to just taking full steps after
  * full_steps_after_k iterations.
  */
class LineSearchFullStepAfterKIter_Step : public LineSearch_Step {
public:

  /// <<std comp>> members for the line search step
  STANDARD_COMPOSITION_MEMBERS(LineSearch_Step,line_search);

  /** \brief . */
  LineSearchFullStepAfterKIter_Step(
        const line_search_ptr_t&	line_search			= 0
      , int						full_steps_after_k
                      = std::numeric_limits<int>::max()	)
    : line_search_(line_search)
      , full_steps_after_k_(full_steps_after_k)
  {}

  /// 
  void full_steps_after_k( int full_steps_after_k )
  {	full_steps_after_k_ = full_steps_after_k; }
  /** \brief . */
  value_type full_steps_after_k() const
  {	return full_steps_after_k_; }

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
  int		full_steps_after_k_;
  
};	// end class LineSearchFullStepAfterKIter_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FULL_STEP_AFTER_K_ITER_STEP_H
