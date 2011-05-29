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

#ifndef COP_QP_SOLVER_STATS_H
#define COP_QP_SOLVER_STATS_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Class for storing statistics about a run of a (active set?) QP solver.
  */
class QPSolverStats {
public:

  // Public types

  /// Set to this value if a statistic is not known.
  enum { NOT_KNOWN = -1 };

  /// Enumeration for the type of point returned from solve_qp(...).
  enum ESolutionType {
    SOLUTION_TYPE_NOT_KNOWN = static_cast<int>(NOT_KNOWN),
    OPTIMAL_SOLUTION		= 0,
    PRIMAL_FEASIBLE_POINT	= 1,
    DUAL_FEASIBLE_POINT		= 2,
    SUBOPTIMAL_POINT		= 3
    };
  /// Enumeration for the type of projected QP on output
  enum EConvexity {
    CONVEXITY_NOT_KNOWN = static_cast<int>(NOT_KNOWN),
    CONVEX              = 0,
    NONCONVEX           = 1
  };

  // Public interface

  /// Construct all unknowns
  QPSolverStats()
    : solution_type_(SOLUTION_TYPE_NOT_KNOWN)
    , convexity_(CONVEXITY_NOT_KNOWN)
    , num_qp_iter_(NOT_KNOWN)
    , num_adds_(NOT_KNOWN), num_drops_(NOT_KNOWN)
    , warm_start_(false), infeasible_qp_(false)
  {}
  /// Initialize the statistics
  void set_stats(
    ESolutionType solution_type, EConvexity convexity
    ,int num_qp_iter, int num_adds, int num_drops
    , bool warm_start, bool infeasible_qp )
  {
    solution_type_	= solution_type;
    convexity_      = convexity;
    num_qp_iter_	= num_qp_iter; 
    num_adds_		= num_adds;
    num_drops_		= num_drops;
    warm_start_		= warm_start;
    infeasible_qp_	= infeasible_qp;
  }
  /** \brief . */
  ESolutionType solution_type() const
  {
    return solution_type_;
  }
  /** \brief . */
  EConvexity convexity() const
  {
    return convexity_;
  }
  /** \brief . */
  int num_qp_iter() const
  {
    return num_qp_iter_;
  }
  /** \brief . */
  int	num_adds() const
  {
    return num_adds_;
  }
  /** \brief . */
  int	num_drop() const
  {
    return num_drops_;
  }
  /** \brief . */
  int	warm_start() const
  {
    return warm_start_;
  }
  /** \brief . */
  int	infeasible_qp() const
  {
    return infeasible_qp_;
  }

private:
  ESolutionType	solution_type_;
  EConvexity      convexity_;
  int				num_qp_iter_;
  int				num_adds_;
  int				num_drops_;
  bool			warm_start_;
  bool			infeasible_qp_;

};	// end class QPSolverStats

inline
std::string toString( const QPSolverStats::ESolutionType &solution_type )
{
  switch(solution_type) {
    case QPSolverStats::SOLUTION_TYPE_NOT_KNOWN:
      return "SOLUTION_TYPE_NOT_KNOWN";
      break;
    case QPSolverStats::OPTIMAL_SOLUTION:
      return "OPTIMAL_SOLUTION";
      break;
    case QPSolverStats::PRIMAL_FEASIBLE_POINT:
      return "PRIMAL_FEASIBLE_POINT";
      break;
    case QPSolverStats::DUAL_FEASIBLE_POINT:
      return "DUAL_FEASIBLE_POINT";
      break;
    case QPSolverStats::SUBOPTIMAL_POINT:
      return "SUBOPTIMAL_POINT";
      break;
    default:
      TEST_FOR_EXCEPT(true);
  }
}

}	// end namespace ConstrainedOptPack

#endif	// COP_QP_SOLVER_STATS_H
