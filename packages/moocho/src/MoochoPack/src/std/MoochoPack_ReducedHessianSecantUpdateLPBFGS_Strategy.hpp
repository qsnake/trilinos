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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_H

#include "MoochoPack_ReducedHessianSecantUpdateBFGSProjected_Strategy.hpp"
#include "MoochoPack_BFGSUpdate_Strategy.hpp"
#include "MoochoPack_quasi_newton_stats.hpp"
#include "MoochoPack_act_set_stats.hpp"
#include "ConstrainedOptPack_MatrixHessianSuperBasic.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Perform BFGS updates on only the free independent (super basic) variables.
 *
 * This method should be very efficient for few super basic variables.
 */
class ReducedHessianSecantUpdateLPBFGS_Strategy : public ReducedHessianSecantUpdate_Strategy
{
public:
  
  /** \brief <<std comp>> members for the strategy object that will
   * perform dense projected BFGS updating.
   */
  STANDARD_COMPOSITION_MEMBERS( ReducedHessianSecantUpdateBFGSProjected_Strategy, proj_bfgs_updater );

  /** \brief Set the minimum number of BFGS updates to perform on the LBFGS matrix
   * before considering switching to projected BFGS updating.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, min_num_updates_proj_start );

  /** \brief Set the maximum number of BFGS updates to perform on the LBFGS matrix
   * before automatically switching to the projected BFGS updating
   * reguardless if the active set has calmed down or not.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, max_num_updates_proj_start );

  /** \brief Set the maximum number of superbasic variables under which switching
   * from limited memory to dense projected PBFGS updating will be allowed.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_superbasics_switch_dense );

  /** \brief Set maximum number of previous BFGS updates to initialize the new dense
   * projected BFGS matrix with.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_add_recent_updates );

  /** \brief . */
    ReducedHessianSecantUpdateLPBFGS_Strategy(
    const proj_bfgs_updater_ptr_t&  proj_bfgs_updater            = NULL
    ,size_type                      min_num_updates_proj_start   = 0
    ,size_type                      max_num_updates_proj_start   = 999999
    ,size_type                      num_superbasics_switch_dense = 500
    ,size_type                      num_add_recent_updates       = 10
    );

  /** \brief . */
  bool perform_update(
    DVectorSlice* s_bfgs, DVectorSlice* y_bfgs, bool first_update
    ,std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,MatrixOp *rHL_k
    );
  /** \brief . */
  void print_step( std::ostream& out, const std::string& leading_str ) const;

private:
  
  // //////////////////////////////
  // Private types

  // /////////////////////////////
  // Private data members

  quasi_newton_stats_iq_member    quasi_newton_stats_;
  act_set_stats_iq_member         act_set_stats_;

}; // end class ReducedHessianSecantUpdateLPBFGS_Strategy

}  // end namespace MoochoPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_LPBFGS_STRATEGY_H
