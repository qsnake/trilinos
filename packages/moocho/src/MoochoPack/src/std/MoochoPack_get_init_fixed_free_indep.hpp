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

#ifndef GET_INIT_FIXED_FREE_INDEP_H
#define GET_INIT_FIXED_FREE_INDEP_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Determine the set of initially fixed and free independent variables.
 *
 * This function will drop all but one of the fixed independent variables
 * whos Lagrange multiplier is above a predefined value.
 *
 * ToDo: Finish documentation.
 *
 * @param  n         [in] Total number of variables.
 * @param  r         [in] Number of decomposed constraints.
 * @param  nu_indep  [in] Sparse vector (size == n-r) of Lagrange multipliers for
 *                   the independent variables.
 * @param  super_basic_mult_drop_tol
 *                   [in] Tolerance of nu_indep(i)/||nu_indep||inf below which
 *                   active variables will not be dropped from the superbasis.
 * @param  olevel    [in] Printing level.
 * @param  out       [out] Stream the output is printed to based on olevel.
 * @param  n_pz_X    [out] Number of dropped super basic variables (n_pz_X <= nu_indep.nz()).
 * @param  n_pz_R    [out] Number of free super basic variables (n-r == n_pz_R + n_pz_X)
 * @param  i_x_free  [out] Array (size >= n_pz_R) of indices the free independent (superbasic) variables.
 * @param  i_x_fixed [out] Array (size >= n_pz_X) of indices the droped (nonbasic) variables.
 * @param  bnd_fixed [out] Array (size >= n_pz_X) of the bounds for the dropped (nonbasic) variables.
 */
void get_init_fixed_free_indep(
  const size_type                        n
  ,const size_type                       r
  ,const SpVectorSlice                   &nu_indep
  ,const value_type                      super_basic_mult_drop_tol
  ,EJournalOutputLevel                   olevel
  ,std::ostream                          &out
  ,size_type                             *n_pz_X
  ,size_type                             *n_pz_R
  ,size_type                             i_x_free[]
  ,size_type                             i_x_fixed[]
  ,ConstrainedOptPack::EBounds  bnd_fixed[]
  );

} // end namespace MoochoPack

#endif // GET_INIT_FIXED_FREE_INDEP_H
