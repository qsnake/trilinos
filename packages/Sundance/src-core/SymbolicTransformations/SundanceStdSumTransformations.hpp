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

#ifndef SUNDANCE_STDSUMTRANSFORMATION_H
#define SUNDANCE_STDSUMTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "SundanceSumTransformationSequence.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;




/**
 * Apply a standard set of transformations
 */
class StdSumTransformations : public SumTransformationSequence
{
public:
  StdSumTransformations();

  virtual ~StdSumTransformations(){;}
};

/** 
 * Rewrite a sum as a polynomial object, if possible.
 */
class IdentifyPolynomialSum : public SumTransformation
{
public:
  /** */
  IdentifyPolynomialSum() : SumTransformation() {;}

  /** */
  virtual ~IdentifyPolynomialSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Put terms in a standard order
 */
class ReorderSum : public SumTransformation
{
public:
  /** */
  ReorderSum() : SumTransformation() {;}

  /** */
  virtual ~ReorderSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};
    

/** 
 * Simplify sums involving unary minuses
 * \f[
 * x + (-y) \rightarrow x-y. 
 * \f]
 * \f[
 * (-x) + (-y) \rightarrow -(x+y). 
 * \f]
 * \f[
 * (-x) + y \rightarrow y-x. 
 * \f]
 */
class RemoveUnaryMinusFromSum : public SumTransformation
{
public:
  /** */
  RemoveUnaryMinusFromSum() : SumTransformation() {;}

  /** */
  virtual ~RemoveUnaryMinusFromSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};
    
/** 
 * Transform a sum by removing a zero term: 
 * \f[
 * x + 0 \rightarrow x. 
 * \f]
 */
class RemoveZeroFromSum : public SumTransformation
{
public:
  /** */
  RemoveZeroFromSum() : SumTransformation() {;}

  /** */
  virtual ~RemoveZeroFromSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Sum two constant exprs without transformation 
 */
class SumConstants : public SumTransformation
{
public:
  /** */
  SumConstants() : SumTransformation() {;}

  /** */
  virtual ~SumConstants(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Transform a sum by moving any constants to the left:
 * \f[
 * x + a \rightarrow a + x
 * \f]
 **/
class MoveConstantsToLeftOfSum : public SumTransformation
{
public:
  /** */
  MoveConstantsToLeftOfSum() : SumTransformation() {;}

  /** */
  virtual ~MoveConstantsToLeftOfSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};


/** 
 * Rearrange a sum whose right operand is a sum including a constant
 * such that constants are grouped on the left:
 * \f[
 * \alpha + s_1 (\beta + s_2 u) \rightarrow (\alpha + s_1 \beta) + s_1 s_2 u
 * \f]
 * \f[
 * u + s_1 (\alpha + s_2 v) \rightarrow s_1 \alpha + (u + s_1 s_2 v)
 * \f]
 **/
class RearrangeRightSumWithConstant : public SumTransformation
{
public:
  /** */
  RearrangeRightSumWithConstant() : SumTransformation() {;}

  /** */
  virtual ~RearrangeRightSumWithConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Rearrange a sum whose left operand is a sum including a constant
 * such that constants are grouped on the left:
 * \f[
 * (\alpha + s_1 u) + s_2 \beta \rightarrow (\alpha + s_2 \beta) + s_1 u 
 * \f]
 * \f[
 * (\alpha + s_1 u) + s_2 v \rightarrow \alpha + (s_1 u + s_2 v)
 * \f]
 **/
class RearrangeLeftSumWithConstant : public SumTransformation
{
public:
  /** */
  RearrangeLeftSumWithConstant() : SumTransformation() {;}

  /** */
  virtual ~RearrangeLeftSumWithConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};


/** 
 * Transform sum of integrals to SumOfIntegral objects
 **/
class SumIntegrals : public SumTransformation
{
public:
  /** */
  SumIntegrals() : SumTransformation() {;}

  /** */
  virtual ~SumIntegrals(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};
}


#endif
