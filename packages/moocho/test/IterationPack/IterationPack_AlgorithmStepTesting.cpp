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

#include "IterationPack_AlgorithmStepTesting.hpp"
#include "IterationPack_Algorithm.hpp"
#include "IterationPack_print_algorithm_step.hpp"

namespace {
char step_type_name[3][15] = { "DO_MAIN_STEP", "DO_PRE_STEP" , "DO_POST_STEP" };
} // namespace

namespace IterationPack {

bool AlgorithmStepTesting::do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
    , poss_type assoc_step_poss)
{
  print_algorithm_step( algo, step_poss, type, assoc_step_poss
    , algo.track().journal_out() );
  return true;
}

void AlgorithmStepTesting::initialize_step(
  Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  print_step_poss( algo, step_poss, type, assoc_step_poss );
  algo.track().journal_out() << "\" : initialize_step(...) called\n";
}

void AlgorithmStepTesting::inform_updated(
  Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  print_step_poss( algo, step_poss, type, assoc_step_poss );
  algo.track().journal_out() << "\" : inform_step(...) called\n";
}

void AlgorithmStepTesting::finalize_step(
  Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  print_step_poss( algo, step_poss, type, assoc_step_poss );
  algo.track().journal_out() << "\" : finalize_step(...) called\n";
}

void AlgorithmStepTesting::print_step( const Algorithm& algo, poss_type step_poss, EDoStepType type
  , poss_type assoc_step_poss ,std::ostream& out, const std::string& leading_str ) const
{
  algo.track().journal_out()
    << std::endl << leading_str << step_poss << ", " << step_type_name[type];
  if(type == DO_MAIN_STEP) {
    algo.track().journal_out()
      << ", \"" << algo.get_step_name(step_poss) << "\"";
  }
  else {
    EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
    algo.track().journal_out()
      << ", " << assoc_step_poss
      << ", \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
  }
  algo.track().journal_out()
    << "\" : print_step(algo,step_poss,type,assoc_step_poss,out) called\n";
}

// private

void AlgorithmStepTesting::print_step_poss(
  const Algorithm& algo, poss_type step_poss, EDoStepType type
  ,poss_type assoc_step_poss
  ) const
{
  algo.track().journal_out()
    << std::endl << step_poss << ", " << step_type_name[type];
  if(type == DO_MAIN_STEP) {
    algo.track().journal_out()
      << ", \"" << algo.get_step_name(step_poss) << "\"";
  }
  else {
    EAssocStepType _type = (type == DO_PRE_STEP ? PRE_STEP : POST_STEP );
    algo.track().journal_out()
      << ", " << assoc_step_poss
      << ", \"" << algo.get_assoc_step_name(step_poss,_type,assoc_step_poss) << "\"";
  }
}

}	// end namespace IterationPack 
