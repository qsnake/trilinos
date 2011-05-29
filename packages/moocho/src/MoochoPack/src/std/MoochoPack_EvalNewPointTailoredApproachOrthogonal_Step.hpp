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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H

#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implements "orthogonal" decompostion for "Tailored Appraoch".
 *
 * Computes:
 \verbatim
 py = inv(I + D*D') * py
 Y  = [ I; -D' ]
 Uy = ???
 \endverbatim
 */
class EvalNewPointTailoredApproachOrthogonal_Step
  : public EvalNewPointTailoredApproach_Step
{
public:

  /** \brief . */
  EvalNewPointTailoredApproachOrthogonal_Step(
    const deriv_tester_ptr_t                &deriv_tester
    ,const bounds_tester_ptr_t              &bounds_tester
    ,EFDDerivTesting                        fd_deriv_testing = FD_DEFAULT
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

  // ///////////////////////////////
  // Private types

  /** \brief . */
  typedef Teuchos::RCP<MatrixSymOpNonsing>  S_ptr_t;

  // ///////////////////////////////
  // Private data members

  S_ptr_t   S_ptr_;

  // //////////////////////////////
  // Private member functions

  // not defined and not to be called
  EvalNewPointTailoredApproachOrthogonal_Step();

};	// end class EvalNewPointTailoredApproachOrthogonal_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
