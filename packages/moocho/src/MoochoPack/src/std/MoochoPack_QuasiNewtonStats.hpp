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

#ifndef QUASI_NEWTON_STATS_H
#define QUASI_NEWTON_STATS_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Class for storing statistics about the Quasi-Newton updating
  */
class QuasiNewtonStats {
public:

  // Public types

  /// Set to this value if a statistic is not known.
  enum EUpdate { UNKNOWN, REINITIALIZED, UPDATED, DAMPENED_UPDATED
    , SKIPED, INDEF_SKIPED };

  // Public interface

  /// Construct all unknowns
  QuasiNewtonStats()
    : update_(UNKNOWN)
  {}

  /// Initialize the statistics
  void set_updated_stats( EUpdate update )
  {
    update_ = update;
  }

  /** \brief . */
  EUpdate updated() const
  {
    return update_;
  }

private:
  EUpdate update_;

};	// end class QuasiNewtonStats

}	// end namespace MoochoPack

#endif	// QUASI_NEWTON_STATS_H
