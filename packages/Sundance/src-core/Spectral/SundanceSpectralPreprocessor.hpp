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

#ifndef SUNDANCE_SPECTRAL_PREPROCESSOR_H
#define SUNDANCE_SPECTRAL_PREPROCESSOR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{

using namespace Sundance;
using namespace Teuchos;

class ProductExpr;
class SumExpr;
class UnaryMinus;
class DiffOp;
class MultiIndex;
  



/** */
class SpectralPreprocessor
{
public:

  /** */
  static Expr projectSpectral(const Expr& e);

  /** */
  static Expr projectSpectral(const Array<Array<Expr> >& terms);

  /** */
  static void parseProduct(const Array<Expr>& factors,
    Expr& test, Expr& deterministic, Array<Expr>& spectral);

  /** */
  static bool isSpectralTest(const Expr& f);

  /** */
  static bool isSpectral(const Expr& f);

  /** */
  static void expandSpectral(const Expr& e, 
    Array<Array<Expr> >& terms);

  /** */
  static bool hasSpectral(const Expr& e);

  /** */
  static void expandSpectralProduct(const ProductExpr* prod,
    Array<Array<Expr> >& terms);

  /** */
  static void expandSpectralSum(const SumExpr* sum,
    Array<Array<Expr> >& terms);

  /** */
  static void expandSpectralUnaryMinus(const UnaryMinus* u,
    Array<Array<Expr> >& terms);

  /** */
  static void expandSpectralDiffOp(const DiffOp* d,
    Array<Array<Expr> >& terms);

  /** */
  static void expandDerivative(const MultiIndex& mi,
    const Array<Expr>& factors,
    Array<Array<Expr> >& productRuleTerms);

  /** */
  static Expr takeDeriv(const Expr& f, const MultiIndex& mi);

  /** */
  static bool isZeroExpr(const Expr& f);
};
  

}

#endif
