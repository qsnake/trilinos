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

#ifndef PRE_PROCESS_BARRIER_LINE_SEARCH_STEP_H
#define PRE_PROCESS_BARRIER_LINE_SEARCH_STEP_H

#include <list>

#include "MoochoPack_Types.hpp"
#include "MoochoPack_LineSearchFilter_Step.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Fraction to boundary rule for calculating alpha max
 *   
 *
 * This class updates alpha_vl, alpha_vu, and alpha
 *  and x_kp1, Vl_kp1, and Vu_kp1
 *
 */

class PreProcessBarrierLineSearch_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
  {
  public:

    /** @name Constructors / initializers */
    //@{

    /** \brief Fraction to Boundary parameter
     *
     * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tau_boundary_frac );

    /** \brief Constructor.
     */
    PreProcessBarrierLineSearch_Step(
      Teuchos::RCP<NLPInterfacePack::NLPBarrier> barrier_nlp,
      const value_type tau_boundary_frac = 0.99
      );
    //@}

    /** @name Overridden from AlgorithmStep */
    //@{
    /** \brief . */
    bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
           , poss_type assoc_step_poss);
    
    
    void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
             , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
    //@}

  private:

    // //////////////////////////
    // Private data members

    Teuchos::RCP<NLPInterfacePack::NLPBarrier> barrier_nlp_;
    CastIQMember< Filter_T > filter_;

  }; // end class PreProcessBarrierLineSearch_Step


class PreProcessBarrierLineSearch_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode,
    public OptionsFromStreamPack::SetOptionsToTargetBase< PreProcessBarrierLineSearch_Step >
  {
  public:
    PreProcessBarrierLineSearch_StepSetOptions(
      PreProcessBarrierLineSearch_Step* target = 0,
      const char opt_grp_name[] = "PreProcessBarrierLineSearch" );

  protected:

    /// Overridden from SetOptionsFromStreamNode
    void setOption( int option_num, const std::string& option_value );
  
  };	// end class PreProcessBarrierLineSearch_StepSetOptions

} // end namespace MoochoPack

#endif // #if !defined PRE_PROCESS_BARRIER_LINE_SEARCH_STEP_H
