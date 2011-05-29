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

#ifndef ALGORITHM_TRACK_COMPOSITE_H
#define ALGORITHM_TRACK_COMPOSITE_H

#include <list>

#include "IterationPack_AlgorithmTracker.hpp"
#include "Teuchos_RCP.hpp"

namespace IterationPack {

/** \brief This class acts a composite container for other \c AlgorithmTracker objects.
 *
 * This class exposes a <tt>std::list<AlgorithmTracker*></tt> object and lets the client
 * manipulate the list.  It is up to the client to maintain this list.
 *
 * See the "Composite" pattern in "Design Patterns", Gama et al, 1995.
 */
class AlgorithmTrackerComposite : public AlgorithmTracker {
public:

  /** \brief . */
  typedef Teuchos::RCP<AlgorithmTracker>      track_ptr_t;
  /** \brief . */
  typedef std::list<track_ptr_t>                                    track_list_t;
  /** \brief . */
  AlgorithmTrackerComposite(const ostream_ptr_t& journal_out);
  /// Give access to the list of \c AlgorithmTracker object pointers.
  track_list_t& tracks();
  /** \brief . */
  const track_list_t& tracks() const;

  /**  @name Overridden from AlgorithmTracker */
  //@{

  /** \brief . */
  void initialize();
  /** \brief . */
  void output_iteration(const Algorithm& algo) const;
  /** \brief . */
  void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  AlgorithmTracker  *tracks;
#else
  track_list_t    tracks_;
#endif

};	// end class AlgorithmTrackerComposite

// ///////////////////////////////////
// Inline members

inline
AlgorithmTrackerComposite::track_list_t&
AlgorithmTrackerComposite::tracks()
{ 
  return tracks_;
}

inline
const AlgorithmTrackerComposite::track_list_t&
AlgorithmTrackerComposite::tracks() const
{ 
  return tracks_;
}

}	// end namespace IterationPack 

#endif	// ALGORITHM_TRACK_COMPOSITE_H
