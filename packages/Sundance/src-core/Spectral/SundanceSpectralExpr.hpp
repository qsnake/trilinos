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

#ifndef SUNDANCE_SPECTRALEXPR_H
#define SUNDANCE_SPECTRALEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceFuncSetAccumulator.hpp"

namespace Sundance
{
using namespace Teuchos;
  
  
  

/**
 * Spectral Expression 
 */

class SpectralExpr : public FuncSetAccumulator, public ScalarExpr
{
private:
  Array<Expr> coeffs_;
  SpectralBasis sbasis_;

public:
  /** Constructor*/
  SpectralExpr (const SpectralBasis& sbasis, const Array<Expr>& coeffs);
  /** Constructor*/
  SpectralExpr (const SpectralBasis& sbasis, const Expr& coeffs);
    
  /** virtual destructor */
  virtual ~SpectralExpr() {;}



  /** Return the Spectral Basis */
  SpectralBasis getSpectralBasis() const;


  /** Return the coefficient of the nth basis term */
  Expr getCoeff(int n) const;

  /** */
  Expr spectralDotProduct(const SpectralExpr* other) const ;

  /** */
  virtual XMLObject toXML() const ;

  /** Write self in text form */
  virtual std::ostream& toText(std::ostream& os, bool paren) const;
 
  /** */
  void accumulateFuncSet(Set<int>& funcDofIDs, 
    const Set<int>& activeSet) const ;

  /** */
  bool hasTestFunctions() const ;

  /** */
  bool hasUnkFunctions() const ;

  /** */
  bool hasHungryDiffOp() const ;

  /** Ordering operator for use in transforming exprs 
   * to standard form */
  bool lessThan(const ScalarExpr* other) const ;


  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
    

};

/** \relates Expr */
Expr getSpectralCoeff(int i, const Expr& e);

/** */
SpectralBasis getSpectralBasis(const Expr& e) ;

}

#endif
