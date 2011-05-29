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

#ifndef PRE_EVAL_NEW_POINT_BARRIER_STEP_H
#define PRE_EVAL_NEW_POINT_BARRIER_STEP_H


#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Standard evaluation step class for extra parameters in primal/dual barrier method.
 *
 * This class calculates \c invXu, \c invXl \c invXu_m_invXl
 *
 */

class PreEvalNewPointBarrier_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
  {
  public:

    /** \brief relative fraction for initializing x within
     *   bounds.
     *   xl_sb = min(xl+relative_bound_push*(xu-xl),
     *               xl + absolute_bound_push)
     *   xu_sb = max(xu-relative_bound_push*(xu-xl),
     *               xu - absolute_bound_push)
     *   if (xl_sb > xu_sb) then
     *      x = (xl + (xu-xl)/2
     *   else if (x < xl_sb) then 
     *      x = xl_sb
     *   else if (x > xu_sb) then
     *      x = xu_sb
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, relative_bound_push );

    /** \brief absolute measure for initializing x within
     *   bounds.
     *   xl_sb = min(xl+relative_bound_push*(xu-xl),
     *               xl + absolute_bound_push)
     *   xu_sb = max(xu-relative_bound_push*(xu-xl),
     *               xu - absolute_bound_push)
     *   if (xl_sb > xu_sb) then
     *      x = (xl + (xu-xl)/2
     *   else if (x < xl_sb) then 
     *      x = xl_sb
     *   else if (x > xu_sb) then
     *      x = xu_sb
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, absolute_bound_push );

    /** @name Overridden from AlgorithmStep */
    //@{
    /** \brief . */
    bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
           , poss_type assoc_step_poss);
    
    
    void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
             , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
    //@}

    /** Constructor.
     */
    PreEvalNewPointBarrier_Step(
      const value_type relative_bound_push = 0.01,
      const value_type absolute_bound_push = 0.001
      );
    //@}

  }; // end class PreEvalNewPointBarrier_Step

class PreEvalNewPointBarrier_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode,
    public OptionsFromStreamPack::SetOptionsToTargetBase< PreEvalNewPointBarrier_Step >
  {
  public:
    PreEvalNewPointBarrier_StepSetOptions(
      PreEvalNewPointBarrier_Step* target = 0,
      const char opt_grp_name[] = "PreEvalNewPointBarrier" );

  protected:

    /// Overridden from SetOptionsFromStreamNode
    void setOption( int option_num, const std::string& option_value );
  
  };	// end class PreEvalNewPointBarrier_StepSetOptions


}  // end namespace MoochoPack

#endif // PRE_EVAL_NEW_POINT_BARRIER_STEP_H
