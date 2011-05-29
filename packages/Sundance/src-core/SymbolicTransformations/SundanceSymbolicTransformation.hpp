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

#ifndef SUNDANCE_SYMBOLICTRANSFORMATION_H
#define SUNDANCE_SYMBOLICTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
using namespace Sundance;


  

class ScalarExpr;

/** 
 * SumTransformation is a base class for any transformation
 * which takes the two operands of a sum (left, right) and produces
 * a new expression mathematically equivalent to the original
 * sum. This will be used to effect simplification
 * transformations on sum expressions.
 */
class SymbolicTransformation : public ObjectWithClassVerbosity<SymbolicTransformation>
{
public:
  /** */
  SymbolicTransformation();

  /** */
  virtual ~SymbolicTransformation(){;}

  /** Whether to build optimized polynomial objects */
  static bool& useOptimizedPolynomials() 
    {static bool rtn=false; return rtn;}

  /** Returns -expr if sign == -1, otherwise returns expr */
  static RCP<ScalarExpr> chooseSign(int sign, 
    const RCP<ScalarExpr>& expr);

  /** Returns -expr if sign == -1, otherwise returns expr */
  static Expr chooseSign(int sign, 
    const Expr& expr);

  /** extract the underlying ScalarExpr from an Expr. */
  static RCP<ScalarExpr> getScalar(const Expr& expr); 
};
}



#endif
