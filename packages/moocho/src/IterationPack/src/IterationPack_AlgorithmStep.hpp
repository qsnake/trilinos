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

#ifndef ALGORITHM_STEP_H
#define ALGORITHM_STEP_H

#include <iosfwd>

#include "IterationPack_Types.hpp"

namespace IterationPack {

/** \brief Base type for all objects that perform steps in an <tt>Algorithm</tt>.
  */
class AlgorithmStep {
public:

  /** \brief . */
  typedef size_t poss_type;

  /** @name Pure virtual functions that must be overridden */
  //@{

  /** \brief Called by <tt>Algorithm</tt> to perform a main, pre or post step at step_poss and assoc_step_poss.
    *
    * @return Should return false if this step object has terminated the algorithm or
    * redirected control to another step.  In this case it is assumed that <tt>this</tt> called
    * <tt>algo\ref Algorithm::terminate ".terminate(...)"</tt>
    * or <tt>algo\ref Algorithm::do_step_next ".do_step_next(...)"</tt>.
    */
  virtual bool do_step(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    ) = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief . */
  virtual ~AlgorithmStep() {}

  /** \brief Called by <tt>Algorithm</tt> just before the algorithm is run.
   *
   * This allows step objects to reinitialize themselves just before
   * an algorithm is run.
   *
   * The default implementation does nothing.
   */
  virtual void initialize_step(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    )
    {}

  /** \brief Called by <tt>Algorithm</tt> to inform when a runtime configuration change
   * is finihed.
   *
   * This function is only called when the algorithm is already running
   * but the configuration has changed.
   *
   * The default implementation does nothing.
   */
  virtual void inform_updated(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    )
    {}

  /** \brief Called by <tt>Algorithm</tt> just after an algorithm is terminiated.
   *
   * This allows step objects to perform any final processing or cleanup
   * just after an algorithm is finished.
   *
   * The default implementation does nothing.
   */
  virtual void finalize_step(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    )
    {}

  /** \brief Called by <tt>Algorithm::print_algorithm()</tt> to print out what this step does in Matlab like format.
    *
    * The default does nothing.
    */
  virtual void print_step(
    const Algorithm&         algo
    ,poss_type               step_poss
    ,EDoStepType             type
    ,poss_type               assoc_step_poss
    , std::ostream&          out
    ,const std::string&      leading_str
    ) const
  {}

  //@}

};	// end class AlgorithmStep

}	// end namespace IterationPack 

#endif // ALGORITHM_STEP_H
