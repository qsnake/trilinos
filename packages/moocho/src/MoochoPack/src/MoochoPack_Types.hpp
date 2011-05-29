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

#ifndef REDUCED_SPACE_SQP_PACK_TYPES_H
#define REDUCED_SPACE_SQP_PACK_TYPES_H

#include "ConstrainedOptPack_Types.hpp"
#include "IterationPack_Types.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerbosityLevel.hpp"

namespace MoochoPack {

using Teuchos::RCP;

// using types from ConstrainedOptPack
#include "ConstrainedOptPack_PublicTypes.ud"

// using types from IterationPack
#include "IterationPack_PublicTypes.ud"

/** \brief enum for journal output. */
enum EJournalOutputLevel {
  PRINT_NOTHING = 0
  ,PRINT_BASIC_ALGORITHM_INFO = 1
  ,PRINT_ALGORITHM_STEPS = 2
  ,PRINT_ACTIVE_SET = 3
  ,PRINT_VECTORS = 4
  ,PRINT_ITERATION_QUANTITIES = 5
};

/** \brief Conver to Teuchos::EVerbosityLevel. */
inline Teuchos::EVerbosityLevel convertToVerbLevel( const EJournalOutputLevel output_level )
{
  switch(output_level) {
    case PRINT_NOTHING:
    case PRINT_BASIC_ALGORITHM_INFO:
      return Teuchos::VERB_NONE;
    case PRINT_ACTIVE_SET:
    case PRINT_ALGORITHM_STEPS:
      return Teuchos::VERB_LOW;
    case PRINT_VECTORS:
      return Teuchos::VERB_HIGH;
    case PRINT_ITERATION_QUANTITIES:
      return Teuchos::VERB_EXTREME;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  return Teuchos::VERB_NONE; // Should never be called!
}

// public interface classes

class NLPAlgoState;
class NLPSolverClientInterface;
class NLPAlgoClientInterface;
class NLPAlgoConfig;

//

class NLPAlgo;

}	// end namespace MoochoPack 

#endif // REDUCED_SPACE_SQP_PACK_TYPES_H
