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

#ifndef RSQP_ALGO_CLIENT_INTERFACE_H
#define RSQP_ALGO_CLIENT_INTERFACE_H

#include "MoochoPack_NLPSolverClientInterface.hpp"

namespace MoochoPack {

/** \brief Interface that smart clients use to set the algorithm configuration
 * object that defines the rSQP algorithm to be used to solve the NLP.
 *
 * ToDo: Finish documentation!
 */
class NLPAlgoClientInterface : public NLPSolverClientInterface {
public:

  /** @name Public Types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<NLPAlgoConfig>	config_ptr_t;

  //@}

  /** @name «std comp» members for config. */
  //@{

  /** \brief . */
  virtual void set_config(const config_ptr_t& config) = 0;
  /** \brief . */
  virtual config_ptr_t& get_config() = 0;
  /** \brief . */
  virtual const config_ptr_t& get_config() const = 0;
  /** \brief . */
  virtual NLPAlgoConfig& config() = 0;
  /** \brief . */
  virtual const NLPAlgoConfig& config() const = 0;

  //@}
  
  /** \brief Causes the algorithm to be configured.
   *
   * Causes the \c config object to configure the algorithm
   * to be ready to solve an NLP or print the algorithm.
   *
   * May be called after the \c nlp, \c track and \c config objects
   * are set.
   *
   * Must be  called before \c print_algorithm() or \c find_min() are called.
    */
  virtual void configure_algorithm(std::ostream* trase_out = 0) = 0;

  /// Print the configured algorithm
  virtual void print_algorithm(std::ostream& out) const = 0;

private:

#ifdef DOXYGEN_COMPILE // Strictly for doxygen diagrams
  /** \brief . */
  NLPAlgoConfig    *config;
#endif

};	// end class NLPAlgoClientInterface

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_CLIENT_INTERFACE_H
