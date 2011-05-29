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

#ifndef SUNDANCE_TEMPSTACK_H
#define SUNDANCE_TEMPSTACK_H

#include "SundanceDefs.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceNoncopyable.hpp"
#include <stack>

namespace Sundance
{
using namespace Sundance;
/**
 * TempStack provides a stack of temporary variables for use during
 * evaluation. 
 *
 * During the course of evaluating an expression, it is often necessary
 * to create temporary variables. For example, in evaluating
 * \code
 * a += b*(c+d)
 * \endcode
 * it is required to create three temporaries (ignoring for
 * explanatory purposes any copies made in cases where one of
 * the vectors will be used elsewhere). We can see this
 * by breaking the operation
 * down into the following steps:
 * \code
 * 1. Create a temporary variable t1
 * 2. Evaluate expression b into t1
 * 3. Create a temporary variable t2
 * 4. Evaluate expression c into t2
 * 3. Create a temporary variable t3
 * 4. Evaluate expression d into t3
 * 5. Carry out t2 += t3
 * 6. Carry out t1 *= t2
 * 7. Carry out a += t1
 * \endcode
 * The number of temporaries required for a given expression
 * will depend on the graph of the expression. In general, we want to
 * create exactly as many temporaries as are needed, and reuse any
 * temporaries that are no longer needed. This is a well-known problem
 * in compiler design, and can be accomplished by maintaining a
 * stack of temporaries. When a new temporary is needed, it is popped
 * from the stack; if the stack is empty, a new temporary is allocated.
 * When a step of a calculation is done, any temporaries used are
 * put back on the stack for further use.
 */
class TempStack : public Noncopyable
{
public:
  /** Empty ctor */
  TempStack();

  /** Construct with an initial vector size */
  TempStack(int vecSize);

  /** Push vector data onto the stack */
  void pushVectorData(const RCP<Array<double> >& vecData) ;

  /** Pop vector data from the stack */
  RCP<Array<double> > popVectorData() ;

  /** Get a new vector (which will often reuse stack data) */
  RCP<EvalVector> popVector() 
    {return rcp(new EvalVector(this));}

  /** */
  void setVecSize(int vecSize) {vecSize_ = vecSize;}

  /** */
  void resetCounter() ;

  /** */
  int numVecsAccessed() const {return numVecsAccessed_;}

  /** */
  int numVecsAllocated() const {return numVecsAllocated_;}

  /** */
  int vecSize() const {return vecSize_;}

private:
          
  int vecSize_;

  std::stack<RCP<Array<double> > > stack_;

  int numVecsAllocated_;

  int numVecsAccessed_;
};
}

#endif
