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

#include <algorithm>

#include "IterationPack_AlgorithmTrackerComposite.hpp"

namespace IterationPack {

AlgorithmTrackerComposite::AlgorithmTrackerComposite(const ostream_ptr_t& journal_out)
  : AlgorithmTracker(journal_out)
{}

void AlgorithmTrackerComposite::initialize()
{
  track_list_t::const_iterator
    itr = tracks_.begin(), itr_end = tracks_.end();
  for(; itr != itr_end; ++itr)
    (*itr)->initialize();
}

void AlgorithmTrackerComposite::output_iteration(
  const Algorithm& algo
  ) const
{
  track_list_t::const_iterator
    itr = tracks_.begin(), itr_end = tracks_.end();
  for(; itr != itr_end; ++itr)
    (*itr)->output_iteration(algo);
}

void AlgorithmTrackerComposite::output_final(
  const Algorithm& algo, EAlgoReturn algo_return
  ) const
{
  track_list_t::const_iterator
    itr = tracks_.begin(), itr_end = tracks_.end();
  for(; itr != itr_end; ++itr)
    (*itr)->output_final(algo,algo_return);
}

} // end namespace IterationPack
