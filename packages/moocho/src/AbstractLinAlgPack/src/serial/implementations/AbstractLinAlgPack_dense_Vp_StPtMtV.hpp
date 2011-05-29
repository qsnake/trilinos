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

#ifndef DENSE_V_P_S_T_P_T_M_T_V_H
#define DENSE_V_P_S_T_P_T_M_T_V_H

#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "MiWorkspacePack.h"

namespace AbstractLinAlgPack {

/** \brief Implements y = b*y + a*op(P)*op(M)*x where it is assumed that
 * generating rows of op(M) is cheap compared to computing op(M)*x.
 *
 * This function is not ment to be used by clients directly but instead
 * by #MatrixOp# subclasses that want to use it.
 */
template<class M_t, class V_t>
void dense_Vp_StPtMtV(
  DVectorSlice                 *y
  ,value_type                 a
  ,const GenPermMatrixSlice   &P
  ,BLAS_Cpp::Transp           P_trans
  ,const M_t                  &M
  ,BLAS_Cpp::Transp           M_trans
  ,const V_t                  &x
  ,value_type                 b
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  using DenseLinAlgPack::dot;
  using DenseLinAlgPack::DVector;
  using DenseLinAlgPack::Vt_S;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::GenPermMatrixSlice;
  typedef AbstractLinAlgPack::EtaVector eta;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  
  // Validate the sizes
  // 
  // y += op(P)*op(M)*x
  // 
  const DenseLinAlgPack::size_type
    ny = y->size(),
    nx = x.size(),
    opM_rows = rows( M.rows(), M.cols(), M_trans ),
    opM_cols = cols( M.rows(), M.cols(), M_trans );
  if(    ny != rows( P.rows(), P.cols(), P_trans )
       || nx != opM_cols
       || cols( P.rows(), P.cols(), P_trans ) != opM_rows )
    throw std::length_error( "MatrixOp::Vp_StPtMtV(...) : Error, "
      "sizes of arguments does not match up" );
  //
  // Compute y = b*y + a*op(P)*op(M)*x in a resonably efficient manner.  This
  // implementation will assume that M is stored as a dense matrix.  Either
  // t = op(M)*x is computed first (O(opM_rows*nx) flops) then y = b*y + a*op(P)*t
  // (O(ny) + O(P_nz) flops) or each row of t' = e(j)' * op(M) (O(nx) flops)
  // is computed one at a time and then y(i) = b*y(i) + a*t'*x (O(nx)  flops)
  // where op(P)(i,j) = 1.0.  In the latter case, there are P_nz rows
  // of op(M) that have to be generated so the total cost is O(P_nz*nx).
  // Therefore, we will do the former if opM_rows < P_nz and the latter otherwise.
  //
  if( !P.nz() ) {
    // y = b*y
    if(b==0.0)       *y = 0.0;
    else if(b!=1.0)  Vt_S(y,b);
  }
  else if( opM_rows > P.nz() || P.is_identity() ) {
    // t = op(M)*x
    Workspace<value_type> t_ws(wss,opM_rows);
    DVectorSlice t(&t_ws[0],t_ws.size());
    LinAlgOpPack::V_MtV( &t, M, M_trans, x );
    // y = b*y + a*op(P)*t
    Vp_StMtV( y, a, P, P_trans, t(), b );
  }
  else {
    // y = b*y
    if(b==0.0)       *y = 0.0;
    else if(b!=1.0)  Vt_S(y,b);
    // Compute t' = e(j)' * op(M) then y(i) += a*t'*x where op(P)(i,j) = 1.0
    Workspace<value_type> t_ws(wss,opM_cols);
    DVectorSlice t(&t_ws[0],t_ws.size());
    if( P.is_identity() ) {
      for( size_type k = 1; k <= P.nz(); ++k ) {
        const size_type
          i = k,
          j = k;
        // t = op(M') * e(j)			
        LinAlgOpPack::V_MtV( &t, M, trans_not(M_trans), eta(j,opM_rows)() );
        // y(i) += a*t'*x
        (*y)(i) += a * dot( t(), x );
      }
    }
    else {
      for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
        const DenseLinAlgPack::size_type
          i = P_trans == no_trans ? itr->row_i() : itr->col_j(),
          j = P_trans == no_trans ? itr->col_j() : itr->row_i();
        // t = op(M') * e(j)			
        LinAlgOpPack::V_MtV( &t, M, trans_not(M_trans), eta(j,opM_rows)() );
        // y(i) += a*t'*x
        (*y)(i) += a * dot( t(), x );
      }
    }
  }
}

} // end namespace AbstractLinAlgPack

#endif // DENSE_V_P_S_T_P_T_M_T_V_H
