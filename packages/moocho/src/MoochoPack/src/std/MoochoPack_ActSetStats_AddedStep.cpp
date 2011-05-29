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

#include "MoochoPack_ActSetStats_AddedStep.hpp"
#include "MoochoPack_active_set_change.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"

bool MoochoPack::ActSetStats_AddedStep::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  if( s.nu().updated_k(0) ) {
    size_type
      num_active = 0, num_adds = 0, num_drops = 0,
      num_active_indep = 0, num_adds_indep = 0, num_drops_indep = 0;
    const SpVector &nu_k = s.nu().get_k(0);
    num_active = nu_k.nz();
    if( s.nu().updated_k(-1) ) {
      const SpVector &nu_km1	= s.nu().get_k(-1);
      active_set_change(
        nu_k(), nu_km1(), s.var_indep(), olevel, &out
        ,&num_adds, &num_drops, &num_active_indep, &num_adds_indep, &num_drops_indep );
      act_set_stats_(s).set_k(0).set_stats(num_active,num_adds,num_drops
                         ,num_active_indep,num_adds_indep,num_drops_indep);
    }
    else {
      act_set_stats_(s).set_k(0).set_stats(
        num_active, ActSetStats::NOT_KNOWN, ActSetStats::NOT_KNOWN
        , nu_k(s.var_indep()).nz(), ActSetStats::NOT_KNOWN, ActSetStats::NOT_KNOWN );
    }
  }
  else {
    if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out	<< "\nnu not calculated for the kth iteration\n";
    }
  }

  return true;
}

void MoochoPack::ActSetStats_AddedStep::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Updates active set statistics for changes from the last iteration\n"
    << L << "Given nu_km1 and nu_k update:\n"
    << L << "    act_set_stats_k(num_active,num_adds,num_drops,num_active_indep,num_adds_indep,num_drops_indep)\n";
}

#endif // 0
