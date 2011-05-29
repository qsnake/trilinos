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

#ifndef SUNDANCE_EVALCONTEXT_H
#define SUNDANCE_EVALCONTEXT_H


#include "SundanceDefs.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Utils.hpp"
#include <algorithm>


namespace Sundance
{
using namespace Teuchos;
using namespace Sundance;
using Sundance::Set;

/** 
 * Different contexts might require the same expression to be
 * evaluated to different orders of functional differentiation; for
 * example, in setting up a linear system, second-order derivatives
 * are required, but in evaluating a functional only zeroth derivs
 * are required. 
 * An EvaluationContext is used as a key to associate an evaluator and
 * its corresponding set of
 * functional derivatives with a context.
 *
 * They key consists of three parts: first, an integer identifier
 * indicating the caller, e.g., an assembler or functional evaluator,
 * second, a set indicating which orders of 
 * differentiation are required by the top level caller, and third,
 a region-quadrature combination.  
*/
class EvalContext
{
public:
  /** Empty ctor */
  EvalContext() : data_() {;}

  /** Construct with a region-quadrature combination and
   * an identifier of the construcing context. */
  EvalContext(const RegionQuadCombo& rqc,
    const Set<int>& needsDiffOrder,
    int contextID)
    : setupVerbosity_(0), 
      maxDiffOrder_(*std::max_element(needsDiffOrder.begin(), needsDiffOrder.end())),
      data_(rcp(new OrderedTriple<Set<int>, int, RegionQuadCombo>(needsDiffOrder, contextID, rqc)))
    {}

  /** Set the verbosity level to be used during preprocessing 
   * of expressions in this context */
  void setSetupVerbosity(int v) {setupVerbosity_ = v;}

  /** Return the verbosity level to be used during preprocessing 
   * of expressions in this context */
  int setupVerbosity() const {return setupVerbosity_;}

  /** Comparison operator for use in maps */
  bool operator<(const EvalContext& other) const 
    {return *data_ < *other.data_;}
          
  /** Write to a std::string */
  std::string toString() const
    {return "EvalContext[diffOrder=" 
        + Teuchos::toString(data_->a())
        + ", id=" 
        + Teuchos::toString(data_->b())
        + ", " + data_->c().toString() + "]";}
          
  /** Write a short description to a std::string */
  std::string brief() const
    {return "EvalContext[diffOrder=" 
        + Teuchos::toString(data_->a())
        + ", id=" 
        + Teuchos::toString(data_->b())
        + "]";}

  /** */
  int topLevelDiffOrder() const {return maxDiffOrder_;}

  /** Indicate whether or not a given order of differentiation 
   * is needed in this context */
  bool needsDerivOrder(int order) const {return data_->a().contains(order);}
  

  /** Return a unique context ID */
  static int nextID() {static int rtn=0; return rtn++;}
private:
  int setupVerbosity_;
  int maxDiffOrder_;
  RCP<OrderedTriple<Set<int>, int, RegionQuadCombo> > data_;
};

}


namespace std
{
/** \relates Sundance::EvalContext */
inline ostream& operator<<(std::ostream& os, 
  const Sundance::EvalContext& c)
{
  os << c.toString();
  return os;
}
}

namespace Teuchos
{
/** \relates Sundance::EvalContext */
inline std::string toString(const Sundance::EvalContext& h)
{return h.toString();}

}


#endif
