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

#include "MoochoPack_ReducedHessianSecantUpdateBFGSFull_Strategy.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

ReducedHessianSecantUpdateBFGSFull_Strategy::ReducedHessianSecantUpdateBFGSFull_Strategy(
  const bfgs_update_ptr_t&      bfgs_update
  )
  :bfgs_update_(bfgs_update)
{}

bool ReducedHessianSecantUpdateBFGSFull_Strategy::perform_update(
  VectorMutable           *s_bfgs
  ,VectorMutable          *y_bfgs
  ,bool                   first_update
  ,std::ostream           & out
  ,EJournalOutputLevel    olevel
  ,NLPAlgo                *algo
  ,NLPAlgoState           *s
  ,MatrixSymOp            *rHL_k
  )
{
  bfgs_update().perform_update(
    s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
    ,rHL_k, &quasi_newton_stats_(*s).set_k(0)
    );
  return true;
}

void ReducedHessianSecantUpdateBFGSFull_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Perform BFGS update on full matrix where: B = rHL_k\n";
  bfgs_update().print_step(out,L);
}

}  // end namespace MoochoPack
