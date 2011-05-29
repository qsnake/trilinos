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

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"

namespace ConstrainedOptPack {

/** \brief Implementation of initial KKT system for all variables initially free
 * and <tt>Ko = G</tt>.
 *
 * In this implementation, \c G must support the \c MatrixSymOpNonsing
 * interface.  Using this initial KKT system essentially make QPSchur use
 * a range space approach (Nocedal & Wright, 1999) for factorizing the KKT
 * system for the current active set.
 */
class QPSchurInitKKTSystemHessianFull
  : public QPSolverRelaxedQPSchur::InitKKTSystem 
{
public:

  /** @name Overridden from InitKKTSystem */
  //@{

  /** \brief Initialize the KKT system where all variables (except the relaxation variable)
   * are initially free and no constraints are in Ko.
   *
   * For this implementation:
   *
   * <tt>n_R = nd</tt>\\
   * <tt>i_x_free = emply (it is identity)</tt>\\
   * <tt>i_x_fixed[0] = nd+1</tt>\\
   * <tt>bnd_fixed[0] = LOWER</tt>\\
   * <tt>j_f_decomp[] = empty</tt>\\
   * <tt>b_X = etaL</tt>\\
   * <tt>Ko = G</tt>\\
   * <tt>fo = -g</tt>\\
   */
  void initialize_kkt_system(
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
    ) const;

  //@}

}; // end class QPSchurInitKKTSystemHessianFull

} // end namesapce ConstrainedOptPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
