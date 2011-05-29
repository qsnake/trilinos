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

#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT

#ifndef QP_SOLVER_RELAXED_QPOPT_H
#define QP_SOLVER_RELAXED_QPOPT_H

#include "ConstrainedOptPack_QPSolverRelaxedQPOPTSOL.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief QPSolver subclass that uses QPOPT.
 *
 * ToDo: Finish documentation.
 */  
class QPSolverRelaxedQPOPT : public QPSolverRelaxedQPOPTSOL
{
public:

  /** \brief . */
  typedef QPSolverRelaxedQPOPTSOL inherited;

  /** \brief Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac );

  /** \brief . */
  QPSolverRelaxedQPOPT(
    value_type max_qp_iter_frac	= 10.0
    );

  /** \brief . */
  ~QPSolverRelaxedQPOPT();

  // /////////////////////////////////
  // Overridden from QPSolverRelaxed

  /** \brief . */
  void release_memory();

protected:

  // /////////////////////////////////////////////////////////////
  // Overridden from QPSolverRelaxedQPOPTSOL

  /** \brief . */
  f_int liwork(f_int N, f_int NCLIN) const;
  /** \brief . */
  f_int lrwork(f_int N, f_int NCLIN) const;
  /** \brief . */
  EInform call_qp_solver(bool warm_start);

private:

  // ////////////////////////////
  // Private types

  /** \brief . */
  enum EQPOPTInform {
    STRONG_LOCAL_MIN      = 0,
    WEAK_LOCAL_MIN        = 1,
    UNBOUNDED             = 2,
    INFEASIBLE            = 3,
    ITMAX_EXCEEDED        = 4,
    MAX_DOF_TOO_SMALL     = 5,
    INVALID_INPUT         = 6,
    PROB_TYPE_NOT_REGOG   = 7
  };

  // ////////////////////////////
  // Private data members

  // extra QPOPT control and input parameters.

  // control

  f_int        ITMAX_;
  f_dbl_prec   BIGBND_;
  f_dbl_prec   FEATOL_;

  // input/output

  f_int           LDA_;
  f_int           LDH_;
  f_dbl_prec*     H_;
  f_int           INFORM_;

  // ////////////////////////////
  // Private member functions

  // not defined and not to be called.
  QPSolverRelaxedQPOPT(const QPSolverRelaxedQPOPT&);
  QPSolverRelaxedQPOPT& operator=(const QPSolverRelaxedQPOPT&);

};	// end class QPSolverRelaxedQPOPT

}	// end namespace ConstrainedOptPack

#endif // QP_SOLVER_RELAXED_QPOPT_H

#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
