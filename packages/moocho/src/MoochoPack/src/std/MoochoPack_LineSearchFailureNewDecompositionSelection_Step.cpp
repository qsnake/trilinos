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
#include <typeinfo>

#include "MoochoPack_LineSearchFailureNewDecompositionSelection_Step.hpp"
#include "MoochoPack_MoochoAlgorithmStepNames.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "IterationPack_print_algorithm_step.hpp"

namespace MoochoPack {

LineSearchFailureNewDecompositionSelection_Step::LineSearchFailureNewDecompositionSelection_Step(
  const line_search_step_ptr_t        &line_search_step
  ,const new_decomp_strategy_ptr_t    &new_decomp_strategy
  )
  :line_search_step_(line_search_step)
  ,new_decomp_strategy_(new_decomp_strategy)
  ,last_ls_failure_k_(-100) // has not failed yet
{}

bool LineSearchFailureNewDecompositionSelection_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  )
{
  try {
    return line_search_step().do_step(_algo,step_poss,type,assoc_step_poss);
  }
  catch(const LineSearchFailure& lsf_excpt) {
    NLPAlgo        &algo = rsqp_algo(_algo);
    NLPAlgoState   &s    = algo.rsqp_state();

    EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
    std::ostream& out = algo.track().journal_out();

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
      out	<< "\nLine search failed "
        << " (k = " << algo.state().k() << ")\n"
        << "LineSearchFailure description: " << lsf_excpt.what() << "\n";
    }

    if( last_ls_failure_k_ == s.k() - 1 ) {
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
        out	<< "\nThe line search failed again even with a new decomposition!"
          << " (k = " << algo.state().k() << ")\n"
          << "We quit!\n";
      }
      throw;
    }

    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
      out	<< "\nSelecting a new decomposition..."
        << " (k = " << algo.state().k() << ")\n";
    }

    last_ls_failure_k_ = s.k();
    return new_decomp_strategy().new_decomposition(algo,step_poss,type,assoc_step_poss);
  }
  return false;	// will never be executed.
}

void LineSearchFailureNewDecompositionSelection_Step::print_step(
  const Algorithm& algo
  ,poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  ,std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "do line search step : " << typeName(line_search_step()) << std::endl;
  line_search_step().print_step(algo,step_poss,type,assoc_step_poss,out,L+"    ");
  out
    << L << "end line search step\n"
    << L << "if thrown line_search_failure then\n"
    << L << "  if line search failed at the last iteration also then\n"
    << L << "    throw line_search_failure\n"
    << L << "  end\n"
    << L << "  new decomposition selection : " << typeName(new_decomp_strategy()) << std::endl
    ;
  new_decomp_strategy().print_new_decomposition(
    rsqp_algo(algo),step_poss,type,assoc_step_poss,out, L + "    " );
  out
    << L << "  end new decomposition selection\n"
    << L << "end\n"
    ;
}

} // end namespace MoochoPack
