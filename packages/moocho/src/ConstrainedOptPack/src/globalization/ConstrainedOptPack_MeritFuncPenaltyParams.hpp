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

#ifndef MERIT_FUNC_PENALTY_PARAMS_H
#define MERIT_FUNC_PENALTY_PARAMS_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace ConstrainedOptPack {

/** \brief This class provides interface for setting and retrieving a penalty parameter
  * that many merit functions use {abstract}.
  */
class MeritFuncPenaltyParams {
public:

  /** \brief . */
  class CanNotResize : public std::logic_error
  {public: CanNotResize(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  virtual ~MeritFuncPenaltyParams() {}

  /** @name To be overridden by subclasses */
  //@{

  /** \brief Set the vector space for \c to use for the penalty parameters.
    */
  virtual void set_space_c( const VectorSpace::space_ptr_t& space_c ) = 0;

  /// Get the vector of penalty parameters for setting them
  virtual VectorMutable& set_mu() = 0;

  /// Get the vector of penalty parameters for viewing them
  virtual const Vector& get_mu() const = 0;

  //@}

};	// end class MeritFuncPenaltyParams

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_PENALTY_PARAMS_H
