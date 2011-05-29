/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_FUNCSETACCUMULATOR_H
#define SUNDANCE_FUNCSETACCUMULATOR_H


#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"

namespace Sundance
{
using Sundance::Set;

/**
 * This abstract class defines the interface for accumulating 
 * sets of functions that are available in an expression.
 */
class FuncSetAccumulator 
{
public:

  /** Append to the set of func IDs present in this expression. 
   * \param funcDofIDs the active functions available from this expression
   * \param activeSet the functions required by the client
   * Most classes will not contain functions, so the 
   * default implementation is a no-op. This will be overridden by 
   * function subtypes
   */
  virtual void accumulateFuncSet(Set<int>& funcDofIDs, 
    const Set<int>& activeSet) const {;}

private:
};

}


#endif
