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

#include "IterationPack_print_algorithm_step.hpp"

void IterationPack::print_algorithm_step( const Algorithm& algo
  , Algorithm::poss_type step_poss, EDoStepType type
  , Algorithm::poss_type assoc_step_poss, std::ostream& out )
{
  out << "\n(" << algo.state().k() << ") " << step_poss;
  if(type == DO_MAIN_STEP) {
    out
      << ": \"" << algo.get_step_name(step_poss) << "\"";
  }
  else {
    EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
    int num_assoc_steps = algo.num_assoc_steps(step_poss,_type);
    out << ".";
    switch(_type) {
      case PRE_STEP:
        out << - num_assoc_steps + ((int)assoc_step_poss - 1);
        break;
      case POST_STEP:
        out << assoc_step_poss;
        break;
    }
    out	<< ": \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
  }
  out << "\n";
}
