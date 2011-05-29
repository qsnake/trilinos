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

#include <math.h>

#include <ostream>

#include "MoochoPack_CheckSkipBFGSUpdateStd_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"

namespace MoochoPack {

CheckSkipBFGSUpdateStd_Step::CheckSkipBFGSUpdateStd_Step(
  value_type	skip_bfgs_prop_const
  )
  :skip_bfgs_prop_const_(skip_bfgs_prop_const)
{}

bool CheckSkipBFGSUpdateStd_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  , poss_type assoc_step_poss
  )
{
  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  IterQuantityAccess<MatrixSymOp>
    &rHL = s.rHL();
  if( rHL.updated_k(-1) )	{
    bool skip_update = true;
    if( !s.Ypy().updated_k(-1) || !s.Zpz().updated_k(-1)
      || !s.rGL().updated_k(-1) || !s.c().updated_k(-1) )
    {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out
          << "\n*** Warning, rHL_km1 is updated but one or more of the quantities Ypy_km1, Zpz_km1\n"
            ",rGL_km1 and c_km1 are not updated.  Therefore, there is not sufficient information\n"
            "to determine if to skip the BFGS update or not.  Check storage requirements for\n"
            " the above quantities\n"
            "The BFGS wupdate will be skipped in this case ...\n";
      }	
      skip_update = true;
    }
    else {
      // The information exists for the update so determine
      // if we are in the region to perform the BFGS update.
      //
      // Check if we are to skip the update for this iteration
      const value_type
        nrm_rGL_km1 = s.rGL().get_k(-1).norm_2(),
        nrm_c_km1	= s.c().get_k(-1).norm_2(),
        nrm_Zpz_km1	= s.Zpz().get_k(-1).norm_2(),
        nrm_Ypy_km1	= s.Ypy().get_k(-1).norm_2();
      // ratio = (skip_bfgs_prop_const / sqrt(||rGL_km1|| + ||c_km1||)) * ( ||Zpz_km1|| / ||Ypy_km1|| )
      value_type
        ratio = ( skip_bfgs_prop_const() / ::sqrt( nrm_rGL_km1 + nrm_c_km1 ) )
          * ( nrm_Zpz_km1 / nrm_Ypy_km1 );
      // If ratio < 1.0 then skip the update
      skip_update = ratio < 1.0;
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out
          << std::endl
          << "ratio = (skip_bfgs_prop_const/sqrt(||rGL_km1||+||c_km1||))*(||Zpz_km1||/||Ypy_km1||)\n"
          << "      = (" << skip_bfgs_prop_const() << "/sqrt("<<nrm_rGL_km1<<"+"<<nrm_c_km1<<"))\n"
          << "        * ("<<nrm_Zpz_km1<<"/"<<nrm_Ypy_km1<<")\n"
          << "      = " << ratio << std::endl
          << "ratio " << (skip_update ? '<' : '>' ) << " 1\n"
          << (skip_update
            ? "Skipping BFGS update ...\n"
            : "Perform BFGS update if you can ...\n"  );
      }
    }
    if(	skip_update ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out
          << "\nrHL_k = rHL_km1\n";
      }
      const MatrixSymOp &rHL_km1 = rHL.get_k(-1);
      rHL.set_k(0) = rHL_km1;
      quasi_newton_stats_(s).set_k(0).set_updated_stats(
        QuasiNewtonStats::SKIPED );
      if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
        out << "\nrHL_k =\n" << rHL.get_k(0);
      }
    }
  }
  
  return true;
}

void CheckSkipBFGSUpdateStd_Step::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Check if we should do the BFGS update\n"
    << L << "if rHL_km1 is update then\n"
    << L << "    If Ypy_km1, Zpz_km1, rGL_km1, or c_km1 is not updated then\n"
    << L << "        *** Warning, insufficient information to determine if we should\n"
    << L << "        *** perform the update.  Check for sufficient backward storage.\n"
    << L << "        rHL_k = rHL_km1\n"
    << L << "    else\n"
    << L << "        *** Check if we are in the proper region\n"
    << L << "        ratio = (skip_bfgs_prop_const/sqrt(norm(rGL_km1,2)+norm(c_km1,2)))\n"
    << L << "                 * (norm(Zpz_km1,2)/norm(Ypy_km1,2) )\n"
    << L << "        if ratio < 1 then \n"
    << L << "            rHL_k = rHL_km1\n"
    << L << "        end\n"
    << L << "    end\n"
    << L << "end\n";
}

} // end namespace MoochoPack
