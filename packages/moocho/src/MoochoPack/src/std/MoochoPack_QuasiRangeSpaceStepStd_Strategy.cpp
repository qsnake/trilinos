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

#include "MoochoPack_QuasiRangeSpaceStepStd_Strategy.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace MoochoPack {

bool QuasiRangeSpaceStepStd_Strategy::solve_quasi_range_space_step(
  std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
  ,const Vector& xo, const Vector& c_xo, VectorMutable* v
    )
{
  using LinAlgOpPack::V_InvMtV;
  using LinAlgOpPack::V_StMtV;
  const MatrixOpNonsing
    &R_k = s->R().get_k(0);
  VectorSpace::vec_mut_ptr_t
    vy = R_k.space_rows().create_member();
  // vy = inv(R_k) * c_xo
  V_InvMtV( vy.get(), R_k, BLAS_Cpp::no_trans, c_xo );
  // v = -Y_k*vy
  V_StMtV( v, -1.0, s->Y().get_k(0), BLAS_Cpp::no_trans, *vy );

  return true;
}

void QuasiRangeSpaceStepStd_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out << L << "*** Compute the approximate range space step:\n"
    << L << "vy = inv(R_k) * c_xo\n"
    << L << "v = -Y_k*vy\n";
}

} // end namespace MoochoPack
