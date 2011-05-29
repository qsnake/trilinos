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

#ifndef RSQP_TRACK_CONSOLE_STD_H
#define RSQP_TRACK_CONSOLE_STD_H

#include "MoochoPack_quasi_newton_stats.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "StopWatchPack_stopwatch.hpp"

namespace MoochoPack {

/** \brief This rSQP iteration class provides a tablular output suitable for
  * an 80 char wide console.
  * 
  * Specifically, these object produces a table with the following fields:
  * <ul>
  * <li> k : The SQP iteration counter starting at zero generally.
  * <li> f : (iteration quantity f_k)
  * 			The value of the objective function at the current iteration.
  * 			This value may be scaled and the scaling factor will be
  * 			printed out before the table is produced.
  * <li> ||c||s : (iteration quantity feas_kkt_err_k)
  * 			The scaled value of the constraints norm (usually the infinity
  * 			norm) at the current iteration.  This is the value compared
  * 			to the convergence criteria feas_tol (see the step
  * 			class \c CheckConvergenceStd_AddedStep).
  * <li> ||rGL||s : (iteration quantity opt_kkt_err_k)
  * 			The scaled value of the norm (usually the infinity
  * 			norm) of the reduced gradient of the Lagrangian
  * 			at the current iteration.  This is the value compared
  * 			to the convergence criteria opt_tol (see the step
  * 			class \c CheckConvergenceStd_AddedStep).
  * <li> QN : (iteration quantity quasi_newton_stats_k)
  * 			Information about the quasi-Newton update of the reduced
  * 			Hessian rHL_k.
  * 			<ul>
  * 			<li> IN : rHL_k was reinitialized (identity?)
  * 			<li> UP : A standard quasi-Newton update was performed on
  * 				rHL_km1 -> rHL_k
  * 			<li> DU : A dampened quasi-Newton update (BFGS) was performed
  * 				on rHL_km1 -> rHL_k.
  * 			<li> SK : The quasi-Newton update was skipped because the
  * 				current iterate was in the wrong region.
  * 			<li> IS : The quasi-Newton update (BFGS) was skipped because
  * 				it was not positive definite or illdefined.
  * 			</ul>
  *	<li> \#act : (iteration quantity nu_k.nz())
  *				The number of active variable bounds at the current iteration.
  *	<li> ||Ypy||2 : (iteration quantity Ypy_k)
  *				The 2 norm of the Range space (feasibility) contribution to the
  *				full step d.
  *	<li> ||Zpz||2 : (iteration quantity Zpz_k)
  *				The 2 norm of the Null space (optimality) contribution to the
  *				full step d.
  *	<li> ||d||inf : (iteration quantity d_k)
  *				The infinity norm of the step vector for the primal unknowns
  *				x_kp1 = x_k + alpha_k * d_k.
  * </ul>
  * 
  * The above quantities can tell you a lot about the progress of the SQP
  * algorithm.
  * 
  * ToDo: Finish discussion.
  * 
  * After the algorithm is finished the total solution time will be printed as
  * well as if the solution was found or not.  Also, the number of function and
  * gradient evaluations will be printed.  Note that the timer is started from
  * the moment this object is created or when set_output_stream(...) is called.
  */
class MoochoTrackerConsoleStd
  : public IterationPack::AlgorithmTracker
{
public:

  /// Construct with an output stream (console presumably)
  MoochoTrackerConsoleStd(const ostream_ptr_t& o, const ostream_ptr_t& journal_out);

  /// Set the output stream for console outputting.
  void set_output_stream(const ostream_ptr_t& o);

  /// Get the output stream for console outputting.
  const ostream_ptr_t& get_output_stream() const;

  /** @name Overridden from AlgorithmTracker */
  //@{

  /// Restarts the timer
  void initialize();
  /** \brief . */
  void output_iteration(const Algorithm& algo) const;
  /** \brief . */
  void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

  //@}

protected:

  /// Print the top header to the output
  void print_top_header(const NLPAlgoState &s, const NLPAlgo& algo) const;

  /// Print the header to the output
  void print_header(const NLPAlgoState &s, const NLPAlgo& algo) const;

  std::ostream& o() const
  {	return *o_; }

private:

  // ///////////////////////////////////////////
  // Private types

  enum { NUM_PRINT_LINES = 10 };
  
  // ///////////////////////////////////////////
  // Private data members

  ostream_ptr_t                       o_;
  mutable StopWatchPack::stopwatch    timer_;
  mutable int                         printed_lines_;
  quasi_newton_stats_iq_member        quasi_newton_stats_;

  // Static formating info.
  static int		w_i2_;
  static char		ul_i2_[];
  static int		w_i4_;
  static char		ul_i4_[];
  static int		p2_;
  static int		w_p2_;
  static char		ul_p2_[];
  static int		p3_;
  static int		w_p3_;
  static char		ul_p3_[];

  // ///////////////////////////////////////////
  // Private member funcitons

  // Not defined and not to be called
  MoochoTrackerConsoleStd();
};	// end class MoochoTrackerConsoleStd

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_CONSOLE_STD_H
