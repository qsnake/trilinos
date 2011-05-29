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

#ifndef NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H
#define NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Solves the unconstrained QP subproblem: <tt>min  qp_grad' * pz + (1/2) * pz' * rHL * pz</tt>.
  *
  * The solution to this system is just:<br>
  * <tt>pz = inv(rHL) *qp_grad</tt>.
  *
  * If use_qp_correc is false then:<br>
  *   <tt>qp_grad = rGf</tt>
  * else<br>
  *   <tt>qp_grad = rGf + zeta * ZtHLYpy.<br>
  *
  * Then <tt>Zpz = Z * pz</tt>
  */
class TangentialStepWithoutBounds_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** \brief Set the maximum size for ||pz|| dampening.
   *
   * A value of <tt>max_pz_norm <= 0.0</tt> means not to dampen pz!
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_pz_norm );

  /** \brief Set the number of iterations to dampen pz for.
   *
   * A value of <tt>num_pz_damp_iters <= 0</tt> means not to dampen pz for any
   * iterations!
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, num_pz_damp_iters );

  /** \brief . */
  TangentialStepWithoutBounds_Step();

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

};	// end class TangentialStepWithoutBounds_Step

}	// end namespace MoochoPack 

#endif	// NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H
