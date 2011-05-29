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

#ifndef ALGORITHM_TRACK_H
#define ALGORITHM_TRACK_H

#include <iosfwd>

#include "IterationPack_Types.hpp"
#include "Teuchos_RCP.hpp"

namespace IterationPack {

/** \brief Used to ouput iteration results and other information.
  *
  * This interface can be implemented by outside clients of an iterative
  * algorithm to monitor or "track" the progress of the algorithm.
  *
  * ToDo: Write more documentation!
  */
class AlgorithmTracker {
public:

  /** \brief . */
  virtual ~AlgorithmTracker() {}

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<std::ostream>    ostream_ptr_t;

  //@}

  /** @name Constructors */
  //@{

  /** \brief Construct with an output stream for journal_out.
   *
   * Preconditions:<ul>
   * <li> <tt>journal_out.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>).
   * </ul>
   */
  AlgorithmTracker(const ostream_ptr_t& journal_out);
  
  //@}
  
  /** @name Algorithm iteration state notification */
  //@{
  
  /** \brief Reinitialize the track object right before it is used.
   *
   * The default implementation does nothing.
   */
  virtual void initialize();

  /** \brief Output information about an iteration just completed.
    *
    * The default just does nothing.
    */
  virtual void output_iteration(const Algorithm& algo) const;

  /** \brief Output information about a just completed algorithm.
    *
    * The default just does nothing.
    */
  virtual void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

  //@}

  /** @name Journal file access */
  //@{

  /** \brief Set a smart pointer to the journal file.
    */
  virtual void set_journal_out(const ostream_ptr_t& journal_out);

  /** \brief Get the smart pointer to the journal file.
   */
  const ostream_ptr_t& get_journal_out() const;

  /** \brief Return a reference to a <tt>std::ostream</tt> to be used to output debug information 
    * and the like.
    */
  virtual std::ostream& journal_out() const;

  //@}

private:

#ifndef DOXYGEN_COMPILE
  ostream_ptr_t   journal_out_;
#endif

  // not defined and not to be called
  AlgorithmTracker();
  
};	// end class AlgorithmTracker

}	// end namespace IterationPack 

#endif // ALGORITHM_TRACK_H
