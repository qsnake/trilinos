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

#ifndef SUNDANCE_EXPRFIELDWRAPPER_H
#define SUNDANCE_EXPRFIELDWRAPPER_H


#include "SundanceDefs.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceFieldBase.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Teuchos;

    
/**
 *
 */
class ExprFieldWrapper : public FieldBase
{
public:
  /** */
  ExprFieldWrapper(const Expr& expr) ;

  /** virtual dtor */
  virtual ~ExprFieldWrapper(){;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const ;

  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const ;

  /** */
  virtual int numElems() const {return Expr_size_;}

  /** */
  virtual bool isPointData() const {return isPointData_;}

  /* */
  GET_RCP(FieldBase);
  /**
   * Return the cell filter on which this field is defined 
   */
  virtual const CellFilter& domain() const 
    { // here we return only the first element ()
      return discreteSpace_.cellFilters(indices_[0][0]);
    }

public:
  Expr expr_;

  // this field should be unique for all the variables
  const DiscreteFunctionData* df_;

  DiscreteSpace discreteSpace_;

  // ---- this field can be expicitly asked -----
  //RCP<DOFMapBase> map_;

  Array< Array<int> > indices_;

  int Expr_size_;

  bool isPointData_;
};
}



#endif
