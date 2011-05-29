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

#ifndef ACT_SET_STATS_H
#define ACT_SET_STATS_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Class for storing statistics about the changes in the active set
  * of an SQP algorithm
  */
class ActSetStats {
public:

  // Public types

  /// Set to this value if a statistic is not known.
  enum { NOT_KNOWN = -1 };

  // Public interface

  /// Construct all unknowns
  ActSetStats()
    : num_active_(NOT_KNOWN), num_adds_(NOT_KNOWN), num_drops_(NOT_KNOWN)
    , num_active_indep_(NOT_KNOWN), num_adds_indep_(NOT_KNOWN), num_drops_indep_(NOT_KNOWN)
  {}

  /// Initialize the statistics
  void set_stats(
    int num_active, int num_adds, int num_drops
    ,int num_active_indep, int num_adds_indep, int num_drops_indep
    )
  {
    num_active_        = num_active;
    num_adds_          = num_adds;
    num_drops_         = num_drops;
    num_active_indep_  = num_active_indep;
    num_adds_indep_    = num_adds_indep;
    num_drops_indep_   = num_drops_indep;
  }

  /** \brief . */
  int num_active() const
  {
    return num_active_;
  }
  /** \brief . */
  int	num_adds() const
  {
    return num_adds_;
  }
  /** \brief . */
  int	num_drops() const
  {
    return num_drops_;
  }

  /** \brief . */
  int num_active_indep() const
  {
    return num_active_indep_;
  }
  /** \brief . */
  int	num_adds_indep() const
  {
    return num_adds_indep_;
  }
  /** \brief . */
  int	num_drops_indep() const
  {
    return num_drops_indep_;
  }

private:
  int num_active_;
  int	num_adds_;
  int	num_drops_;
  int num_active_indep_;
  int	num_adds_indep_;
  int	num_drops_indep_;

};	// end class ActSetStats

}	// end namespace MoochoPack

#endif	// ACT_SET_STATS_H
