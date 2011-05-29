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

#ifndef ACTIVE_SET_CHANGE_H
#define ACTIVE_SET_CHANGE_H

#include <iosfwd>

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Calculate the change in the active set and output change
  * if asked to.
  *
  * ToDo: Add more description of the output you get.
  *
  * @param	nu_k	[in] Multipliers for variable bounds for iteration k.
  * @param	num_km1	[in] Multipliers for variable bounds for iteration k-1
  * @param	olevel	[in] Specifies the output level.  We have:\\
  *				PRINT_NOTHING : No output is sent to out\\
  *				PRINT_ALGORITHM_STEPS :
  *					Just the number of additions
  *					and deletions to the active set and the total number
  *					of active constraints is output.\\ 
  *				PRINT_ACTIVE_SET : Enumerates
  *					which variable were added and dropped from the active set.\\
  * @param	num_adds
  *                 [out] Gives the total number of variables fixed at a bound
  *					added to the active set.
  * @param	num_drops
  *                 [out] Gives the total number of variables freed from a 
  *					bound and dropped from the active set.
  * @param	num_adds_indep
  *                 [out] Gives the number of independent variables fixed at a bound
  *					added to the active set.
  * @param	num_drops_indep
  *                 [out] Gives the number of independent variables freed from a 
  *					bound and dropped from the active set.
  * @param	out	[O] Target for output.
  */
void active_set_change(
  const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1, Range1D indep
  ,EJournalOutputLevel olevel, std::ostream* out
  ,size_type* num_adds, size_type* num_drops
  ,size_type* num_active_indep, size_type* num_adds_indep, size_type* num_drops_indep
  ); 

}	// end namespace MoochoPack

#endif	// ACTIVE_SET_CHANGE_H
