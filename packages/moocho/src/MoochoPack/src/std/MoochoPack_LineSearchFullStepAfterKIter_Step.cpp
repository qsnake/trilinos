#if 0

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

#include <ostream>

#include "MoochoPack_LineSearchFullStepAfterKIter_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"

bool MoochoPack::LineSearchFullStepAfterKIter_Step::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  NLP			&nlp	= algo.nlp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();
  out << std::boolalpha;

  // print step header.
  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  const bool take_full_step = s.k() > full_steps_after_k();

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nk = " << s.k() << ( take_full_step ? " > " : " < ")
          << "full_steps_after_k = " << full_steps_after_k() << std::endl;
  }

  if( !take_full_step ) {
    return line_search().do_step(_algo,step_poss,type,assoc_step_poss);
  }
  else {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out	<< "\nKeep the full step...\n";
    }
  }

  return true;
}

void MoochoPack::LineSearchFullStepAfterKIter_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out	<< L << "*** Start using full steps after full_steps_after_k iterations.\n"
    << L << "default: full_steps_after_k = very big\n";
  out	<< L << "if k < full_steps_after_k then\n";
  line_search().print_step(algo,step_poss,type,assoc_step_poss,out,L + "    " );
  out	<< L << "end\n";
}

#endif // 0
