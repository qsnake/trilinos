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

#ifndef MOOCHO_TRACKER_XML_SUMMARY_HPP
#define MOOCHO_TRACKER_XML_SUMMARY_HPP

#include "IterationPack_Algorithm.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "MoochoPack_Types.hpp"

namespace MoochoPack {


using IterationPack::Algorithm;
using IterationPack::EAlgoReturn;

/** \brief This class outputs an XML summary file of the algorithm 
 *   results and performance
 */
class MoochoTrackerXMLSummary
  : public IterationPack::AlgorithmTracker
{
public:

  /// Construct with an output stream
  MoochoTrackerXMLSummary(
    const Teuchos::RCP<std::ostream> &journal_out
    ,const std::string xml_filename
    ,const std::string problem_name
    ,const std::string algorithm_description
    );

  /// Set the output stream for summary outputting
  //void set_output_stream(const ostream_ptr_t& o);

  /// Get the output stream for summary outputting.
  //const ostream_ptr_t& get_output_stream() const;

  /// Output a basic file (with failed status)
  //   that will be overwritten if there is no
  //   exception
  void output_pre_file() const;

  /** @name Overridden from AlgorithmTracker */
  //@{

  /** \brief . */
  void output_iteration(const Algorithm& algo) const;
  /** \brief . */
  void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;
  
  //@}

protected:

  /// Print the header to the output
  void open_problem_element( std::ostream& out, const Algorithm& algo) const;
  void close_problem_element( std::ostream& out) const;

private:

  mutable value_type obj_value_;
  mutable value_type c_norm_value_;

  std::string xml_filename_;
  std::string problem_name_;
  std::string algorithm_description_;

  // Not defined and not to be called
  MoochoTrackerXMLSummary();

};	// end class MoochoTrackerXMLSummary

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_SUMMARY_STD_H
