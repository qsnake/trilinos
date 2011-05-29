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

#ifndef RSQP_ALGO_H
#define RSQP_ALGO_H

#include "MoochoPack_NLPAlgoInterface.hpp"
#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "IterationPack_Algorithm.hpp"
#include "StandardAggregationMacros.hpp"

namespace MoochoPack {

/** \brief rSQP Algorithm control class.
  */
class NLPAlgo
  : public NLPAlgoInterface
  , public IterationPack::Algorithm
{
public:

  /** @name Public Types */
  //@{

  //@}

  /// Constructs with no step, added_step, pre_step, post_step, state, or decomp_sys objects added.
  NLPAlgo();

  /// <<std aggr>> members for algo_cntr
  STANDARD_AGGREGATION_MEMBERS( NLPAlgoContainer, algo_cntr )

  /// <<std aggr>> members for nlp
  STANDARD_AGGREGATION_MEMBERS( NLP, nlp )

  /** \brief . */
  NLPAlgoState& rsqp_state()
  {	return dynamic_cast<NLPAlgoState&>(state()); }

  /** \brief . */
  const NLPAlgoState& rsqp_state() const
  {	return dynamic_cast<const NLPAlgoState&>(state()); }

  /** \brief . */
  void do_step_first(Algorithm::poss_type first_step_poss)
  {	first_step_poss_ = first_step_poss; }

  /** @name Overridden form rSQPAlgoInteface */
  //@{	
  
  /** \brief . */
  const NLPAlgoState& retrieve_state() const;

  /** \brief This is the main control function for the rSQP algorithm.
    *
    * This function basically just calls Algorithm::do_algorithm(...).
    */
  NLPSolverClientInterface::EFindMinReturn dispatch();

  /** \brief . */
  void interface_print_algorithm(std::ostream& out) const;
  /** \brief . */
  void interface_set_algo_timing( bool algo_timing );
  /** \brief . */
  bool interface_algo_timing() const;
  /** \brief . */
  void interface_print_algorithm_times( std::ostream& out ) const;

  //@}

  /// overridden from Algorihth.

  /** \brief . */
  void print_algorithm(std::ostream& out) const;

protected:

  // First step to execute
  Algorithm::poss_type first_step_poss_;

};	// end class NLPAlgo

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_H
