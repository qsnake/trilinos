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

#ifndef EVAL_NEW_POINT_STD_STEP_H
#define EVAL_NEW_POINT_STD_STEP_H

#include "MoochoPack_DecompositionSystemHandler_Strategy.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "ConstrainedOptPack_DecompositionSystemTester.hpp"
#include "ConstrainedOptPack_VariableBoundsTester.hpp"
#include "NLPInterfacePack_NLPFirstDerivTester.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Standard new point evaluation step class.
 *
 * This class calculates \c Gc, \c Gh, updates the range/null decompositon matrices
 * \c Z, \c Y, \c R, \c Uz, \c Uy \c Vz and \c Vy and calculates \c Gf, \c c,
 * \c h, and \c f in that order.
 */
class EvalNewPointStd_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EFDDerivTesting   { FD_DEFAULT,  FD_TEST,  FD_NO_TEST  };

  //@}

  /** @name Constructors / initializers */
  //@{

  /// «std comp» members for range/null decomposition handler
  STANDARD_COMPOSITION_MEMBERS( DecompositionSystemHandler_Strategy, decomp_sys_handler );
  /// «std comp» members for first derivative tester object
  STANDARD_COMPOSITION_MEMBERS( NLPFirstDerivTester, deriv_tester );
  /// «std comp» Members for variable bounds tester object
  STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester );
  /// «std comp» members for decomp_sys tester tester object
  STANDARD_COMPOSITION_MEMBERS( DecompositionSystemTester, decomp_sys_tester );
  /** \brief Set how and if finite derivatives are tested.
    *
    * ToDo: Finish documentation.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDDerivTesting, fd_deriv_testing );
  /** \brief Set how and if the decomposition system is tested.
    *
    * ToDo: Finish documentation.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DecompositionSystemHandler_Strategy::EDecompSysTesting, decomp_sys_testing );
  /** \brief Set how to set the print level for decomp_sys_tester (only if testing).
    *
    * ToDo: Finish documentation.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DecompositionSystemHandler_Strategy::EDecompSysPrintLevel, decomp_sys_testing_print_level );

  /** \brief Constructor.
   *
   * new_point == true by default.
   */
  EvalNewPointStd_Step(
    const decomp_sys_handler_ptr_t                              &decomp_sys_handler
    ,const deriv_tester_ptr_t                                   &deriv_tester
    ,const bounds_tester_ptr_t                                  &bounds_tester
    ,const decomp_sys_tester_ptr_t                              &decomp_sys_tester
    ,EFDDerivTesting                                            fd_deriv_testing               = FD_DEFAULT
    ,DecompositionSystemHandler_Strategy::EDecompSysTesting     decomp_sys_testing             = DecompositionSystemHandler_Strategy::DST_DEFAULT
    ,DecompositionSystemHandler_Strategy::EDecompSysPrintLevel  decomp_sys_testing_print_level = DecompositionSystemHandler_Strategy::DSPL_USE_GLOBAL
    );

  //@}

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

  // Not defined and not to be called
  EvalNewPointStd_Step();

};	// end class EvalNewPointStd_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_STD_STEP_H
