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

#include "MoochoPack_SetDBoundsStd_AddedStep.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"

namespace MoochoPack {

SetDBoundsStd_AddedStep::SetDBoundsStd_AddedStep()
  :dl_iq_(dl_name)
  ,du_iq_(du_name)
  {}

bool SetDBoundsStd_AddedStep::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  NLPAlgo &algo = rsqp_algo(_algo);
  NLPAlgoState &s = algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
  std::ostream &out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  const Range1D
    var_dep = s.var_dep(),
    var_indep = s.var_indep();

  const Vector
    &x_k = s.x().get_k(0),
    &xl  = algo.nlp().xl(),
    &xu  = algo.nlp().xu();
  VectorMutable
    &dl  = dl_iq_(s).set_k(0),
    &du  = du_iq_(s).set_k(0);
  
  // dl = xl - x_k
  LinAlgOpPack::V_VmV( &dl, xl, x_k );	

  // du = xu - x_k
  LinAlgOpPack::V_VmV( &du, xu, x_k );	

  if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    if(var_indep.size()) {
      out << "\ndl(var_indep)_k = \n" << *dl.sub_view(var_indep);
      out << "\ndu(var_indep)_k = \n" << *du.sub_view(var_indep);
    }
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    if(var_dep.size()) {
      out << "\ndl(var_dep)_k = \n" << *dl.sub_view(var_dep);
      out << "\ndu(var_dep)_k = \n" << *du.sub_view(var_dep);
    }
    out << "\ndl_k = \n" << dl;
    out << "\ndu_k = \n" << du;
  }

  return true;
}

void SetDBoundsStd_AddedStep::print_step(
  const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Set the bounds on d\n"
    << L << "d_bounds_k.l = xl - x_k\n"
    << L << "d_bounds_k.u = xu - x_k\n"
    ;
}

} // end namespace MoochoPack
