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

#ifndef RSQP_ALGO_CONTAINER_H
#define RSQP_ALGO_CONTAINER_H

#include "MoochoPack_NLPAlgoClientInterface.hpp"
#include "MoochoPack_NLPAlgoInterface.hpp"
#include "MoochoPack_NLPAlgoConfig.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implementation for NLPAlgo solver.
 *
 * Acts as a container for NLPAlgo.  This class is hidden from clients
 * by not exposing it to them in header files.
 */
class NLPAlgoContainer : public NLPAlgoClientInterface {
public:

  /** @name Constructors / initializers */
  //@{

  /// Members for <<std comp>> of the algorithm object algo.
  STANDARD_COMPOSITION_MEMBERS( NLPAlgoInterface, algo );

  /// Construct a container with no configuration object set.
  NLPAlgoContainer()
  {}

  //@}

  /** @name Overridden from NLPAlgoClientInterface */
  //@{

  /** \brief . */
  void set_config(const config_ptr_t& config);
  /** \brief . */
  config_ptr_t& get_config();
  /** \brief . */
  const config_ptr_t& get_config() const;
  /** \brief . */
  NLPAlgoConfig& config();
  /** \brief . */
  const NLPAlgoConfig& config() const;

  //@}

  /** @name Overridden from NLPSolverClientInterface */
  //@{

  /** \brief . */
  EFindMinReturn find_min();
  /** \brief . */
  void configure_algorithm(std::ostream* trase_out);
  /** \brief . */
  void print_algorithm(std::ostream& out) const;
  /** \brief . */
  void set_algo_timing( bool algo_timing );
  /** \brief . */
  bool algo_timing() const;
  /** \brief . */
  void print_algorithm_times( std::ostream& out ) const;

  //@}

private:

  config_ptr_t			config_;

  // Assert that the object has been set up properly and throw exception if it has not
  void assert_valid_setup() const;

  // Not defined and not to be called
  NLPAlgoContainer(const NLPAlgoContainer&);
  NLPAlgoContainer& operator=(const NLPAlgoContainer&);

};	// end class NLPAlgoContainer

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_CONTAINER_H
