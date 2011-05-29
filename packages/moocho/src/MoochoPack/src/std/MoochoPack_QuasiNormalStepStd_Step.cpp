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

#include "MoochoPack_QuasiNormalStepStd_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace MoochoPack {

bool QuasiNormalStepStd_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  using BLAS_Cpp::no_trans;
  using AbstractLinAlgPack::Vt_S;
  using AbstractLinAlgPack::V_InvMtV;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;

  NLPAlgo         &algo        = rsqp_algo(_algo);
  NLPAlgoState        &s           = algo.rsqp_state();
  const Range1D    equ_decomp   = s.equ_decomp();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // Get iteration quantities
  IterQuantityAccess<VectorMutable>
    &c_iq   = s.c(),
    &py_iq  = s.py(),
    &Ypy_iq = s.Ypy();
  IterQuantityAccess<MatrixOpNonsing>
    &R_iq = s.R();
  IterQuantityAccess<MatrixOp>
    &Y_iq = s.Y();

  // Solve the system py = - inv(R) * c(equ_decomp)
  VectorMutable &py_k = py_iq.set_k(0);
  V_InvMtV( &py_k, R_iq.get_k(0), no_trans, *c_iq.get_k(0).sub_view(equ_decomp) );
  Vt_S( &py_k, -1.0 );

  // Ypy = Y * py
  V_MtV( &Ypy_iq.set_k(0), Y_iq.get_k(0), no_trans, py_k );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\n||py||   = " << py_iq.get_k(0).norm_inf() << std::endl
      << "\n||Ypy||2 = " << Ypy_iq.get_k(0).norm_2()  << std::endl;
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out	<< "\npy_k =\n"  << py_iq.get_k(0);
    out	<< "\nYpy_k =\n" << Ypy_iq.get_k(0);
  }

  return true;
}

void QuasiNormalStepStd_Step::print_step(
  const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Calculate the range space step\n"
    << L << "py_k = - inv(R_k) * c_k(equ_decomp)\n"
    << L << "Ypy_k = Y_k * py_k\n";
}

} // end namespace MoochoPack
