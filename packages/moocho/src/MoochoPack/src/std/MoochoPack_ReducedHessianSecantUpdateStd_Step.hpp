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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H
#define REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H

#include "MoochoPack_ReducedHessianSecantUpdate_Strategy.hpp"
#include "MoochoPack_quasi_newton_stats.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Updates rHL_k using a secant update.
 *
 * This class will initialize rHL = I if it is not already and will handle
 * transitions to new basis selections by reseting rHL = I.  The actually
 * update is delegated to a strategy object (see #secant_update# below).
 *
 * See the printed step documentation for a description.
 */
class ReducedHessianSecantUpdateStd_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** \brief <<std comp>> members for the strategy interface object that will
   * actually perform the secant update.
   */
  STANDARD_COMPOSITION_MEMBERS( ReducedHessianSecantUpdate_Strategy, secant_update );

  /** \brief . */
  ReducedHessianSecantUpdateStd_Step(
    const secant_update_ptr_t& secant_update = Teuchos::null
    );

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

private:
  enum { NO_BASIS_UPDATED_YET = INT_MIN };
  int                             num_basis_;
  int                             iter_k_rHL_init_ident_;
  quasi_newton_stats_iq_member	quasi_newton_stats_;
  
};	// end class ReducedHessianSecantUpdateStd_Step

}	// end namespace MoochoPack 

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_STD_STEP_H
