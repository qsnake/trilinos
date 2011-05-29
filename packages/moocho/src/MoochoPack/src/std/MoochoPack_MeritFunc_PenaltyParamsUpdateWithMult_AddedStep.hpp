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

#ifndef MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Updates a set of penalty parameters for a merit function as:
  * mu(j) = max( mu(j), |lambda_k(j)| ).
  *
  * This class assumes the merit function supports the interfaces
  * MeritFuncPenaltyParams and MeritFuncNLPDirecDeriv.
  */
class MeritFunc_PenaltyParamsUpdateWithMult_AddedStep
  : public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

  /// <<std comp>> members for merit_func
  STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func);

  /** \brief . */
  MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(const merit_func_ptr_t& merit_func
    , value_type small_mu = 1e-6, value_type min_mu_ratio = 1e-8
    , value_type mult_factor = 1e-4 , value_type kkt_near_sol = 1e-1 );

  // ///////////////////////////////
  // Overridden from AlgorithmStep

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss
    , IterationPack::EDoStepType type, poss_type assoc_step_poss
    , std::ostream& out, const std::string& leading_str ) const;

  // //////////////////////////////////////////////////////////////////////
  // Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

  /** \brief . */
  void small_mu( value_type small_mu );
  /** \brief . */
  value_type small_mu() const;

  /** \brief . */
  void min_mu_ratio( value_type min_mu_ratio );
  /** \brief . */
  value_type min_mu_ratio() const;

  /** \brief . */
  void mult_factor( value_type mult_factor );
  /** \brief . */
  value_type mult_factor() const;

  /** \brief . */
  void kkt_near_sol( value_type kkt_near_sol );
  /** \brief . */
  value_type kkt_near_sol() const;

private:
  bool near_solution_;
  value_type small_mu_;
  value_type min_mu_ratio_;
  value_type mult_factor_;
  value_type kkt_near_sol_;
  value_type norm_inf_mu_last_;
  
};	// end class MeritFunc_PenaltyParamsUpdateWithMult_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H
