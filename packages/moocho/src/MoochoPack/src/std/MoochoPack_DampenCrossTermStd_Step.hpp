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

#ifndef DAMPEN_CROSS_TERM_STD_STEP_H
#define DAMPEN_CROSS_TERM_STD_STEP_H

#include "rSQPAlgo_Step.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Compute a dampening term zeta_k for the cross term w_k such that
 * Gf'*Z*pz <= 0.
 *
 * This condition Gf'*Z*pz <= 0 is needed to ensure descent of many
 * merit functions.
 *
 * This implementation ensures that Gf'*Z*pz <= 0 only if
 * there will not be any active constraints when the reduced QP subproblem
 * is solved (nu_k = 0) and there are no undecomposed equality constraints
 * or if there they are linearly dependent (lambda_k(undecomp_con) = 0).
 * 
 * In particular this implementation computes zeta_k such that:
 * 
 * Gf'*Z*pz <= frac_descent * rGf'inv(B)*rGf
 * 
 * where: 0 < frac_descent < 1
 * 
 * To ensure strong descent (and hopefully deal with the cases where
 * nu_k != 0 and lambda_k(undecomp_con) != 0) the parameter frac_descent
 * is set to frac_descent = 0.9 by default.
 * 
 * The basis derivation goes like this:
 * 
 * ToDo: Finish documentation!
 */
class DampenCrossTermStd_Step : public rSQPAlgo_Step {
public:

  /// «std comp» members for frac_descent
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, frac_descent );

  /** \brief . */
  DampenCrossTermStd_Step(const value_type& frac_descent = 0.9);

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class DampenCrossTermStd_Step

}	// end namespace MoochoPack 

#endif	// DAMPEN_CROSS_TERM_STD_STEP_H
