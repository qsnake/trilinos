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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H

#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"

namespace MoochoPack {

/** \brief Implements "coordinate" decompostion for "Tailored Appraoch".
  *
  * Computes:<br>
  * <tt>py = py</tt><br>
  * <tt>Y = [ I; 0 ]</tt><br>
  */
class EvalNewPointTailoredApproachCoordinate_Step
  : public EvalNewPointTailoredApproach_Step
{
public:

  /** \brief . */
  EvalNewPointTailoredApproachCoordinate_Step(
      const deriv_tester_ptr_t& 	deriv_tester
    , const bounds_tester_ptr_t&	bounds_tester
    , EFDDerivTesting				fd_deriv_testing = FD_DEFAULT
    );

protected:

  /** @name Overridden from EvalNewPointTailoredApproach_Step */
  //@{

  /** \brief . */
  void uninitialize_Y_Uy(
    MatrixOp         *Y
    ,MatrixOp        *Uy
    );
  /** \brief . */
  void calc_py_Y_Uy(
    const NLPDirect       &nlp
    ,const D_ptr_t        &D
    ,VectorMutable        *py
    ,MatrixOp             *Y
    ,MatrixOp             *Uy
    ,EJournalOutputLevel  olevel
    ,std::ostream         &out
    );
  /** \brief . */
  void recalc_py(
    const MatrixOp           &D
    ,VectorMutable           *py
    ,EJournalOutputLevel     olevel
    ,std::ostream            &out
    );
  /** \brief . */
  void print_calc_py_Y_Uy(
    std::ostream& out, const std::string& leading_str
    ) const;

  //@}

private:
  // not defined and not to be called
  EvalNewPointTailoredApproachCoordinate_Step();

};	// end class EvalNewPointTailoredApproachCoordinate_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_COORDINATE_STEP_H
