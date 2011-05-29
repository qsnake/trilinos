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

#ifndef MERIT_FUNC_CALC_NLP_H
#define MERIT_FUNC_CALC_NLP_H

#include "ConstrainedOptPack_MeritFuncCalc.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "StandardAggregationMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Adds the ability to compute <tt>phi(f(x),c(x),h(x))</tt> at \c x
  * directly instead of having to compute f, c and h first.
  * This class uses an aggregate NLP to perform the computations of \a f(x)
  * \a c(x) and \a h(x).
  */
class MeritFuncCalcNLP : public MeritFuncCalc {
public:

  /** @name Constructors/initializers */
  //@{

  /// <<std aggr>> stereotype members for phi.
  STANDARD_CONST_AGGREGATION_MEMBERS( MeritFuncNLP, phi )

  /// <<std aggr>> stereotype members for nlp.
  STANDARD_CONST_AGGREGATION_MEMBERS( NLP, nlp )

  /** \brief . */
  MeritFuncCalcNLP( const MeritFuncNLP* phi = 0, const NLP* nlp = 0 );

  //@}

  /** @name Overridden from MeritFuncCalc */
  //@{

  /** \brief Return the value of the merit function at x.
   * Here phi(x) is calculated directly using the nlp.
   */
  value_type operator()(const Vector& x) const;

  /// Calls phi().deriv() on phi.
  value_type deriv() const;

  /// Calls <tt>phi().print_merit_func()</tt>.
  void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

};	// end class MeritFuncCalcNLP

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_NLP_H
