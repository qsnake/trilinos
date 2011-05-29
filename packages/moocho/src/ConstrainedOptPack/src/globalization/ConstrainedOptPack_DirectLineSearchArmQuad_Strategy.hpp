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

#ifndef DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H
#define DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H

#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Performs a line search using the Armijo condition and
  * uses quadratic interpolation to select each new alpha.
  */
class DirectLineSearchArmQuad_Strategy : public DirectLineSearch_Strategy {
public:

  /// Set the Armijo cord test fractional reduction parameter.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta );
  
  /// The minimum fraction that alpha is reduced for each line search iteration.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, min_frac );

  /// The maximum fraction that alpha is reduced for each line search iteration.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_frac );

  /** \brief Deterimine if the line search iterations are maxed out or not.
   * 
   * This option is really only used for debugging and requires
   * changing the other parameters to make it useful.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, max_out_iter );

  /// Constructs with default settings.
  DirectLineSearchArmQuad_Strategy(
    int           max_iter       = 20
    ,value_type   eta            = 1.0e-4
    ,value_type   min_frac       = 0.1
    ,value_type   max_frac       = 0.5
    ,bool         max_out_iter   = false
    );

  /** @name Overridden from DirectLineSearch_Strategy */
  //@{

  /** \brief . */
  void set_max_iter(int max_iter);
  /** \brief . */
  int max_iter() const;
  /** \brief . */
  int num_iterations() const;
  /** \brief Performs the following line search:<br>
   *
   \verbatim

   num_iter = 0;
   while( phi.value(alpha_k) > phi_k + eta * alpha_k * phi.deriv() ) 
   {
      if(num_iter >= max_iter) return true;
      num_iter = num_iter + 1;
      alpha_k = [ min_frac * alpha_k <= quadradic interpolation for alpha	<= max_frac * alpha_k ];
   }
   return true;<br>
   \endverbatim
   * If the maximum number of iterations is exceeded then false will be returned.
   *
   * The default values for the adjustable parameters (from D&S A6.3.1)
   * are:<br>
   * max_iter = 20<br>
   * eta = 1.0e-4<br>
   * min_frac = 0.1<br>
   * max_frac = 0.5<br>
   */
  bool do_line_search(
    const MeritFuncCalc1D   &phi
    ,value_type             phi_k
    ,value_type             *alpha_k
    ,value_type             *phi_kp1
    ,std::ostream           *out
    );

  /** \brief . */
  void print_algorithm(std::ostream& out, const std::string& leading_str) const;

  //@}

private:
  int	max_iter_;
  int	num_iter_;	// stores the number of interations

  // Throw an exception if the parameters are not in a proper range.
  void validate_parameters() const;

};	// end class DirectLineSearchArmQuad_Strategy

}	// end namespace ConstrainedOptPack

#endif	// DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H
