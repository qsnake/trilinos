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

#ifndef SET_D_BOUNDS_STD_ADDED_STEP_HH
#define SET_D_BOUNDS_STD_ADDED_STEP_HH

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"

namespace MoochoPack {

/** \brief Computes the bounds for the QP subproblem from the %NLP bounds.
 */
class SetDBoundsStd_AddedStep
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** \brief . */
  SetDBoundsStd_AddedStep();

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

private:
  IterationPack::CastIQMember<VectorMutable> dl_iq_;
  IterationPack::CastIQMember<VectorMutable> du_iq_;

}; // end class SetDBoundsStd_AddedStep

} // end namespace MoochoPack 

#endif // SET_D_BOUNDS_STD_ADDED_STEP_HH
