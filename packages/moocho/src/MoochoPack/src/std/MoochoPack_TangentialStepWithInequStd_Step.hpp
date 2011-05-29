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

#ifndef NULL_SPACE_STEP_WITH_INEQU_STD_STEP_H
#define NULL_SPACE_STEP_WITH_INEQU_STD_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"
#include "MoochoPack_qp_solver_stats.hpp"
#include "MoochoPack_act_set_stats.hpp"
#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Solves the reduced QP subproblem with bounds and/or general
 * inequalities.
 *
 * ToDo: Finish documentation.
 */
class TangentialStepWithInequStd_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /// QP solver
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxed, qp_solver );

  /// QP solver tester
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxedTester, qp_tester );

  /** \brief Set the ratio of the number of inequality constraints in the
   * active-set of the last two calls before a warm start is attempted.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warm_start_frac );

  /** \brief . */
  enum EQPTesting { QP_TEST_DEFAULT, QP_TEST, QP_NO_TEST };

  /** \brief Set how and if the QP solution is tested.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPTesting, qp_testing );

  /** \brief Determine if a QPFailure exception is thrown if the QP solver
   * returns PRIMAL_FEASIBLE_POINT.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, primal_feasible_point_error );

  /** \brief Determine if a \c QPFailure exception is thrown if the QP solver
   * returns \c DUAl_FEASIBLE_POINT.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dual_feasible_point_error );

  /// Construct and initialize
  TangentialStepWithInequStd_Step(
    const qp_solver_ptr_t       &qp_solver
    ,const qp_tester_ptr_t      &qp_tester
    ,value_type                 warm_start_frac             = 0.8
    ,EQPTesting                 qp_testing                  = QP_TEST_DEFAULT
    ,bool                       primal_feasible_point_error = true
    ,bool                       dual_feasible_point_error   = true
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
  IterationPack::CastIQMember<VectorMutable>  dl_iq_;
  IterationPack::CastIQMember<VectorMutable>  du_iq_;
  qp_solver_stats_iq_member                                qp_solver_stats_;
  act_set_stats_iq_member                                  act_set_stats_;

  // not defined and not to be called
  TangentialStepWithInequStd_Step();

};	// end class TangentialStepWithInequStd_Step

}	// end namespace MoochoPack 

#endif	// NULL_SPACE_STEP_WITH_INEQU_STD_STEP_H
