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

#ifndef SUNDANCE_CELLVECTOREXPR_H
#define SUNDANCE_CELLVECTOREXPR_H

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCellVectorEvaluator.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Sundance;

/** */
enum CellVectorExprType {CellNormalVector, CellTangentSpace};

/**
 * Expression that returns a component of a cell-determined vector for
 * each cell on 
 * which it is evaluated. This makes sense only for boundary cells. 
 */
class CellVectorExpr
  : public EvaluatableExpr,
    public GenericEvaluatorFactory<CellVectorExpr, CellVectorEvaluator>
{
public:
  /** */
  CellVectorExpr(
    int normalComponentIndex, 
    int dimension,
    const std::string& name
    );

  /** */
  CellVectorExpr(
    int tangentBasisIndex, 
    int tangentComponentIndex,
    int dimension,
    const std::string& name
    );


    
  /** */
  virtual ~CellVectorExpr() {;}

  /** */
  bool isTangent() const {return type_ == CellTangentSpace;}

  /** */
  bool isNormal() const {return !isTangent();}

  /** */
  int componentIndex() const {return componentIndex_;}

  /** */
  int basisMemberIndex() const {return basisMemberIndex_;}

  /** */
  int dimension() const {return dim_;}

  /** */
  virtual XMLObject toXML() const ;

  const std::string& name() const {return name_;}

  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;
    

  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;
private:
  std::string name_;
  int dim_;
  CellVectorExprType type_;
  int basisMemberIndex_;
  int componentIndex_;
  
};


/** \relates Expr \relates CellVectorExpr */
Expr CellNormalExpr(int dimension, const std::string& name);

/** \relates Expr \relates CellVectorExpr */
Expr CellTangentExpr(int dimension, const std::string& name);

}

#endif
