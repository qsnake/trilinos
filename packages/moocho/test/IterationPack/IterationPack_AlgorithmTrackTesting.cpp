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

#include <iomanip>

#include "IterationPack_AlgorithmTrackTesting.hpp"
#include "IterationPack_Algorithm.hpp"

namespace IterationPack {

void AlgorithmTrackTesting::output_iteration(const Algorithm& algo) const {
  std::ostream &o = journal_out();
  o
    << "\ntrack.output_iteration(algo) called for the iteration k = "
    << algo.state().k() << std::endl;
  o
    << "\nStep times for the current iteration (in seconds):\n";
  const int n = algo.num_steps();
  std::vector<double> step_times(n+1);
  algo.get_step_times_k(0,&step_times[0]);
  o << "  step_id:time = ";
  for( int l = 0; l < n+1; ++l )
    o << " " << (l+1) << ":" << step_times[l];
  o << "\n  total time = " << step_times[n] << std::endl;
  
}

void AlgorithmTrackTesting::output_final(const Algorithm& algo, EAlgoReturn algo_return) const {
  char algo_return_name[6][50] =
    {
      "TERMINATE_TRUE"
      ,"TERMINATE_FALSE"
      ,"MAX_ITER_EXCEEDED"
      ,"MAX_RUN_TIME_EXCEEDED"
      ,"INTERRUPTED_TERMINATE_TRUE"
      ,"INTERRUPTED_TERMINATE_FALSE"
    };
  std::ostream &o = journal_out();
  o << "\ntrack.output_final(algo,algo_return) called for the iteration k = "
    << algo.state().k() << " and algo_return = " << algo_return_name[algo_return]
    << std::endl;
  o << "Timing (in seconds) statistics for step 0 : ";
  double total, average, min, max, percent;
  algo.get_final_step_stats(0,&total,&average,&min,&max,&percent);
  o << "total = " << total << ", average = " << average << ", min = " << min
    << ", max = " << max << ", percent = " << percent << std::endl;
}

}	// end namespace IterationPack 
