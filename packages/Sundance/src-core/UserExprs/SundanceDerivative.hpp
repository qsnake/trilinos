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

#ifndef SUNDANCE_DERIVATIVE_H
#define SUNDANCE_DERIVATIVE_H

#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceMultiIndex.hpp"
#include "Teuchos_XMLObject.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/** 
 * Derivative is an expression subtype representing 
 * a first-order spatial partial derivative operator.
 */
class Derivative : public ScalarExpr
{
public:
  /** Construct an operator for spatial differentiation with respect to
   * the given direction (0=x, 1=y, or 2=z).  */
  Derivative(int direction);

  /** virtual destructor */
  virtual ~Derivative() {;}

  /** */
  virtual XMLObject toXML() const ;

  /** Indicate whether this expression is a "hungry"
   * differential operator that is awaiting an argument. */
  virtual bool isHungryDiffOp() const {return true;}

  /** */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** */
  const MultiIndex& multiIndex() const {return m_;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


private:
  MultiIndex m_;

};
}

#endif
