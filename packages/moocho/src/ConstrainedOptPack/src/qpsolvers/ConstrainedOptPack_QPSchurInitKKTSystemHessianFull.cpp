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

#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace ConstrainedOptPack {

void QPSchurInitKKTSystemHessianFull::initialize_kkt_system(
  const Vector    &g
  ,const MatrixOp   &G
  ,value_type           etaL
  ,const Vector   *dL
  ,const Vector   *dU
  ,const MatrixOp   *F
  ,BLAS_Cpp::Transp     trans_F
  ,const Vector   *f
  ,const Vector   *d
  ,const Vector   *nu
  ,size_type            *n_R
  ,i_x_free_t           *i_x_free
  ,i_x_fixed_t          *i_x_fixed
  ,bnd_fixed_t          *bnd_fixed
  ,j_f_decomp_t         *j_f_decomp
  ,DVector               *b_X
  ,Ko_ptr_t             *Ko
  ,DVector               *fo
  ) const
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  using LinAlgOpPack::V_mV;

  // Validate type of and convert G
  const MatrixSymOpNonsing&
    G_sym = dyn_cast<const MatrixSymOpNonsing>(G);

  const size_type nd = g.dim();
  
  // n_R
  *n_R = nd;
  // i_x_free[i-1] = i, i = 1...nd
  i_x_free->resize(0);
  // i_x_fixed[0] = nd+1
  i_x_fixed->resize(1);
  (*i_x_fixed)[0] = nd+1;
  // bnd_fixed[0] = LOWER
  bnd_fixed->resize(1);
  (*bnd_fixed)[0] = LOWER;
  // j_f_decomp[] = empty
  j_f_decomp->resize(0);
  // b_X = etaL
  b_X->resize(1);
  (*b_X)[0] = etaL;
  // Ko = G
  *Ko = Teuchos::rcp(&G_sym,false); // Not dynamically allocated so don't delete!
  // fo = -g
  V_mV(fo,VectorDenseEncap(g)());
}

} // end namesapce ConstrainedOptPack
