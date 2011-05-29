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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_FULL_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_FULL_STRATEGY_H

#include "MoochoPack_ReducedHessianSecantUpdate_Strategy.hpp"
#include "MoochoPack_BFGSUpdate_Strategy.hpp"
#include "MoochoPack_quasi_newton_stats.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Perform BFGS updates on full reduced Hessian.
 *
 * This is really a do nothing class that just uses a strategy
 * object (see #bfgs_update# below) to perform the update on the full
 * reduced hessian matrix.
 */
class ReducedHessianSecantUpdateBFGSFull_Strategy : public ReducedHessianSecantUpdate_Strategy
{
public:
  
  /** \brief <<std comp>> members for the strategy object that will
   * perform guts secant update.
   */
  STANDARD_COMPOSITION_MEMBERS( BFGSUpdate_Strategy, bfgs_update );

    ReducedHessianSecantUpdateBFGSFull_Strategy(
    const bfgs_update_ptr_t&      bfgs_update = Teuchos::null
    );      

  /** @name Overridden from ReducedHessianSecantUpdate_Strategy */
  //@{
  /** \brief . */
  bool perform_update(
    VectorMutable     *s_bfgs
    ,VectorMutable    *y_bfgs
    ,bool                   first_update
    ,std::ostream           & out
    ,EJournalOutputLevel    olevel
    ,NLPAlgo               *algo
    ,NLPAlgoState              *s
    ,MatrixSymOp        *rHL_k
    );
  /** \brief . */
  void print_step( std::ostream& out, const std::string& leading_str ) const;
  //@}

private:
  quasi_newton_stats_iq_member	quasi_newton_stats_;

}; // end class ReducedHessianSecantUpdateBFGSFull_Strategy

}  // end namespace MoochoPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_BFGS_FULL_STRATEGY_H
