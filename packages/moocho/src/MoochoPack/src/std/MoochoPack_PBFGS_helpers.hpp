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

#ifndef PBFGS_HELPERS
#define PBFGS_HELPERS

#include "MoochoPack_Types.hpp"
#include "MoochoPack_ActSetStats.hpp"

namespace MoochoPack {
namespace PBFGSPack {

/** @name Helper functions for PBFGS updating.
 *
 */
//@{

/** \brief Determine if the active set has calmed down enough and print this test.
 *
 * This function will return true if:
 \begin{verbatim}
  ( num_adds_indep == NOT_KNOWN || num_drops_indep == NOT_KNOWN || num_active_indep == 0
     ? 0.0
     : std::_MAX(((double)(num_active_indep)-num_adds_indep-num_drops_indep) / num_active_indep, 0.0
   ) >= act_set_frac_proj_start
  &&
  num_active_indep > 0
 \end{verbatim}
 * Otherwise this function will return false.
 */
bool act_set_calmed_down( 
  const ActSetStats         &act_set_stats
  ,const value_type         act_set_frac_proj_start
  ,EJournalOutputLevel      olevel
  ,std::ostream             &out
  );

/** \brief Initialize i_x_free[], s_R'*s_R and s_R'*y_R for free variables not in nu_indep.
 *
 * @param  nu_indep  [in] Sparse vector (n_pz = nu_indep.size(), n_pz_R = n_pz - nu_indep.nz())
 *                   of Lagrange multipliers for the independent variables.
 * @param  s         [in] DVector (size n_pz) secant update vector for B*s=y
 * @param  y         [in] DVector (size n_pz) secant update vector for B*s=y
 * @param  n_pz_R    [out] Number of free super basic variables (n_pz_R = n_pz - nu_indep.nz())
 * @param  i_x_free  [out] Array (size n_pz_R) of indices the free independent (superbasic) variables.
 *                   These are the indices not in nu_indep. 
 * @param  sRTsR     [out] s_R'*s_R
 * @param  yRTyR     [out] y_R'*y_R
 */
void init_i_x_free_sRTsR_sRTyR(
  const SpVectorSlice        &nu_indep
  ,const DVectorSlice         &s
  ,const DVectorSlice         &y
  ,size_type                 *n_pz_R
  ,size_type                 i_x_free[]
  ,value_type                *sRTsR
  ,value_type                *sRTyR
  );

/** \brief Sort fixed variables  according to the condition:
 *
 * #|s_X(i)^2*B(i,i)|/|sRTBRRsR| + |s_X(i)*y_X(i)|/|sRTyR|#.
 *
 * The input quantities are defined as follows:
 \begin{verbatim}
  for k = 0...nu_indep.nz()-1
      i = (nu_indep.begin()+k)->indice() + nu_indep.offset()
      B_ii = B_XX[k]
      s_i = s(i)
      y_i = y(i)
 \end{verbatim}
 *
 * @param  nu_indep  [in] Sparse vector (n_pz = nu_indep.size(), n_pz_R = n_pz - nu_indep.nz())
 *                   of Lagrange multipliers for the independent variables.
 * @param  s         [in] DVector (size n_pz), secant update vector for B*s=y
 * @param  y         [in] DVector (size n_pz), secant update vector for B*s=y
 * @param  B_XX      [in] DVector (size nu_indep.nz()), Diagonal elements for B_XX
 * @param  sRTBRRsR  [in] s_R' * B_RR * s_R
 * @param  sRTyR     [in] s_R' * y_R
 * @param  sXTBXXsX  [out] s_X' * B_XX * s_X
 * @param  sXTyX     [out] s_X' * y_X
 * @param  l_x_fixed_sorted
 *                   [out] Array (size nu_indep.nz()) which gives the indices
 *                   l = l_x_fixed_sorted[k], k = 0...nu_indep.nz()-1, where
 *                   i = (nu_indep.begin() + l)->indice() + nu_indep.offset()
 *                   which are sorted according to: 
 *                   #|s(i)^2*B_XX(l,l)|/|sRTBRRsR| + |s(i)*y(i)|/|sRTyR|#
 */
void sort_fixed_max_cond_viol(
  const SpVectorSlice        &nu_indep
  ,const DVectorSlice         &s
  ,const DVectorSlice         &y
  ,const DVectorSlice         &B_XX
  ,const value_type          sRTBRRsR
  ,const value_type          sRTyR
  ,value_type                *sXTBXXsX
  ,value_type                *sXTyX
  ,size_type                 l_x_fixed_sorted[]
  );

/** \brief Choose the rest of i_x_free[] and i_x_fixed[].
 *
 * The input quantities are defined as follows:
 \begin{verbatim}
  for k = 0...nu_indep.nz()-1
      i = (nu_indep.begin()+k)->indice() + nu_indep.offset()
      B_ii = B_XX[k]
      s_i = s(i)
      y_i = y(i)
 \end{verbatim}
 * This function adjusts i_x_free[] and i_x_fixed[] so that
 * the following conditions are satisfied:
 *
 * #(s_X'*B_XX*s_X)/(s_R'*B_RR*s_R) <= project_error_tol#
 *
 * #(s_X'*y_X')/(s_R'*y_R') <= project_error_tol#
 *
 * #|nu_indep(i_x_fixed[k])|/||nu_indep||inf >= super_basic_mult_drop_tol#
 *
 * @param  project_error_tol
 *                [in] see above
 * @param  super_basic_mult_drop_tol
 *                [in] see above
 * @param  nu_indep  [in] Sparse vector (n_pz = nu_indep.size(), n_pz_R = n_pz - nu_indep.nz())
 *                   of Lagrange multipliers for the independent variables.
 * @param  s         [in] DVector (size n_pz), secant update vector for B*s=y
 * @param  y         [in] DVector (size n_pz), secant update vector for B*s=y
 * @param  B_XX      [in] DVector (size nu_indep.nz()), Diagonal elements for B_XX
 * @param  l_x_fixed_sorted
 *                   [in] Array (size nu_indep.nz()) which gives the indices
 *                   l = l_x_fixed_sorted[k], k = 0...nu_indep.nz()-1, where
 *                   i = (nu_indep.begin() + l)->indice() + nu_indep.offset()
 *                   which are sorted according to:
 *                   #|s(i)^2*B_XX(l,l)|/|sRTBRRsR| + |s(i)*y(i)|/|sRTyR|#
 * @param  olevel    [in] Printing level.
 * @param  out       [out] Stream the output is printed to based on olevel.
 * @param  sRTBRRsR  [in/out] On input, this must equal s_R'*B_RR*s_R where Q_R
 *                   is defined by i_x_free[] on input.  On output, it will equal
 *                   s_R'*B_RR*s_R where Q_R is defined by i_x_free[] on output.
 * @param  sRTyR     [in/out] On input, this must equal s_R'*y_R where Q_R
 *                   is defined by i_x_free[] on input.  On output, it will equal
 *                   s_R'*y_R where Q_R is defined by i_x_free[] on output.
 * @param  sXTBXXsX  [in/out] On input, this must equal s_X'*B_XX*s_X where Q_X
 *                   is defined by what is not in i_x_free[] on input.  On output, it will equal
 *                   s_X_'*B_XX*s_X where Q_X is defined by i_x_fixed[] on output.
 * @param  sXTyX     [in/out] On input, this must equal s_X'*y_X where Q_X
 *                   is defined by what is not in i_x_free[] on input.  On output, it will equal
 *                   s_X_'*y_X where Q_X is defined by i_x_fixed[] on output.
 * @param  n_pz_X    [out] Number of dropped super basic variables (n_pz_X <= nu_indep.nz()).
 * @param  n_pz_R    [in/out] On input contains the number of free variables in i_x_free[] on input.
 *                   On output contaitns the adjusted number of superbasic variables (n-r == n_pz_R + n_pz_X)
 * @param  i_x_free  [in/out] Array (length n_pz_R on output).
 *                   On input i_x_free[0...n_pz_R-1] must contain the indicies as initialized
 *                   by init_i_x_free...(nu_indep,n_pz_R,i_x_free...).  One output it will contain
 *                   any addition indices that define Q_R needed to satsify the above conditions.
 *                   This array must be sorted in accending order on input and will be sorted on output.
 * @param  i_x_fixed [out] Array (length n_pz_X on output).  On output will contain the
 *                   indices that define Q_X that satisfy the above conditions.
 * @param  bnd_fixed [out] Array (lenght n_pz_X on output).  On output will contain the
 *                   bounds of the indices in i_x_fixed[] on output.
 */
void choose_fixed_free(
  const value_type                       project_error_tol
  ,const value_type                      super_basic_mult_drop_tol
  ,const SpVectorSlice                   &nu_indep
  ,const DVectorSlice                     &s
  ,const DVectorSlice                     &y
  ,const DVectorSlice                     &B_XX
  ,const size_type                       l_x_fixed_sorted[]
  ,EJournalOutputLevel                   olevel
  ,std::ostream                          &out
  ,value_type                            *sRTBRRsR
  ,value_type                            *sRTyR
  ,value_type                            *sXTBXXsX
  ,value_type                            *sXTyX
  ,size_type                             *n_pz_X
  ,size_type                             *n_pz_R
  ,size_type                             i_x_free[]
  ,size_type                             i_x_fixed[]
  ,ConstrainedOptPack::EBounds  bnd_fixed[]
  );

//@}

} // end namespace PBFGSPack
} // end namespace MoochoPack

#endif // PBFGS_HELPERS
