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

// #define RELEASE_TRACE

#include <iostream> // for debugging Release version.
#include <typeinfo>

#include "MoochoPack_NLPAlgo.hpp"

namespace MoochoPack {

NLPAlgo::NLPAlgo()
  : algo_cntr_(NULL), nlp_(NULL), first_step_poss_(1)
{}

// Overridden form rSQPAlgoInteface

const NLPAlgoState& NLPAlgo::retrieve_state() const
{
  return dynamic_cast<const NLPAlgoState&>(state());
}

NLPSolverClientInterface::EFindMinReturn
NLPAlgo::dispatch() {
  switch( do_algorithm(first_step_poss_) ) {
    case IterationPack::TERMINATE_TRUE:
      return NLPSolverClientInterface::SOLUTION_FOUND;
    case IterationPack::TERMINATE_FALSE:
      return NLPSolverClientInterface::ALGORITHMIC_ERROR;
    case IterationPack::MAX_ITER_EXCEEDED:
      return NLPSolverClientInterface::MAX_ITER_EXCEEDED;
    case IterationPack::MAX_RUN_TIME_EXCEEDED:
      return NLPSolverClientInterface::MAX_RUN_TIME_EXCEEDED;
    case IterationPack::INTERRUPTED_TERMINATE_TRUE:
      return NLPSolverClientInterface::SOLUTION_FOUND;
    case IterationPack::INTERRUPTED_TERMINATE_FALSE:
      return NLPSolverClientInterface::ALGORITHMIC_ERROR;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return NLPSolverClientInterface::SOLUTION_FOUND;	// will never be called.
}

void NLPAlgo::interface_print_algorithm(std::ostream& out) const {
  print_steps(out);
  print_algorithm(out);
}

void NLPAlgo::interface_set_algo_timing( bool algo_timing ) {
  set_algo_timing(algo_timing);
}

bool NLPAlgo::interface_algo_timing() const {
  return algo_timing();
}

void NLPAlgo::interface_print_algorithm_times( std::ostream& out ) const {
  print_algorithm_times(out);
}

// Overridden from Algorithm.

void NLPAlgo::print_algorithm(std::ostream& out) const {
  out
    << "\n*** NLP ***\n"
    << typeName(*get_nlp()) << "\n";

  Algorithm::print_algorithm(out);
}

}	// end namespace MoochoPack
