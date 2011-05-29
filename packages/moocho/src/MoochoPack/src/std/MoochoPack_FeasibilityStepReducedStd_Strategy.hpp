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

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H

#include "MoochoPack_FeasibilityStep_Strategy.hpp"
#include "MoochoPack_QuasiRangeSpaceStep_Strategy.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implements the feasibility step computation for reduced space SQP.
 */
class FeasibilityStepReducedStd_Strategy : public FeasibilityStep_Strategy
{
public:

  /// <<std comp>> members for the qp solver
  STANDARD_COMPOSITION_MEMBERS( QuasiRangeSpaceStep_Strategy, quasi_range_space_step );

  typedef ConstrainedOptPack::QPSolverRelaxedTester
    QPSolverRelaxedTester;

  /// QP solver
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxed, qp_solver );

  /// Comparision object compatible with Gc
  STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxedTester, qp_tester );
    
  /** \brief . */
  enum EQPObjective {
    OBJ_MIN_FULL_STEP           ///< min 1/2 * (Y*wy + Z*wz)'*(Y*wy + Z*wz)
    ,OBJ_MIN_NULL_SPACE_STEP    ///< min 1/2 * wz'*wz
    ,OBJ_RSQP                   ///< min qp_grad_k'*wz + 1/2 * wz'*rHL_k*wz
  };

  /** \brief Set what is used for the QP objective.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPObjective, qp_objective );

  /** \brief . */
  enum EQPTesting {
    QP_TEST_DEFAULT     ///< Decide based on olevel input to <tt>compute_feasibility_step(...)</tt>
    ,QP_TEST            ///< Perform the tests
    ,QP_NO_TEST         ///< Don't perform the tests
  };

  /** \brief Set how and if the QP solution is tested.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPTesting, qp_testing );

  /// Construct and initialize
  FeasibilityStepReducedStd_Strategy(
    const quasi_range_space_step_ptr_t   &quasi_range_space_step
    ,const qp_solver_ptr_t               &qp_solver
    ,const qp_tester_ptr_t               &qp_tester
    ,EQPObjective                        qp_objective     = OBJ_MIN_NULL_SPACE_STEP
    ,EQPTesting                          qp_testing       = QP_TEST_DEFAULT
    );

  // ////////////////////////////////////////////
  // Overridden from FeasibilityStep_Strategy

  /** \brief Computes a feasibility step by computing simple quasi-range and null space components.
   *
   * ToDo: Finish documentation!
   *
   */
   bool compute_feasibility_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const Vector& xo, const Vector& c_xo, VectorMutable* w
      );

  /** \brief . */
  void print_step( std::ostream& out, const std::string& leading_str ) const;

private:

  IterationPack::CastIQMember<VectorMutable>  dl_iq_;
  IterationPack::CastIQMember<VectorMutable>  du_iq_;
  int                                                      current_k_;
  Teuchos::RCP<const MatrixOp>            Hess_ptr_;
  VectorSpace::vec_mut_ptr_t                               grad_store_;
  DMatrix                                                Hess_store_;

}; // end class FeasibilityStepReducedStd_Strategy

} // end namespace MoochoPack

#endif // FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
