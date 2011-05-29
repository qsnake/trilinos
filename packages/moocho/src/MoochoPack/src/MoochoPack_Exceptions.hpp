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

#ifndef REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
#define REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H

#include "MoochoPack_Types.hpp"
#include "ConstrainedOptPack_QPSolverStats.hpp"

namespace MoochoPack {

/** \defgroup MoochoPack_grp Standard exceptions for MoochoPack */
//@{

// Thrown if the constraints are infeasible
class InfeasibleConstraints : public std::logic_error
{public: InfeasibleConstraints(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a line search failure occurs.
class LineSearchFailure : public std::runtime_error
{public: LineSearchFailure(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a runtime test failed.
class TestFailed : public std::runtime_error
{public: TestFailed(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a the QP failed and was not corrected
class QPFailure : public std::runtime_error
{
public:
  QPFailure(const std::string& what_arg
        , const ConstrainedOptPack::QPSolverStats& _qp_stats)
    : std::runtime_error(what_arg)
    , qp_stats(_qp_stats)
    {}
  ConstrainedOptPack::QPSolverStats qp_stats;
};

//@}

}	// end namespace MoochoPack 

#endif // REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
