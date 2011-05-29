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

#ifndef SUNDANCE_TESTEVALMEDIATOR_H
#define SUNDANCE_TESTEVALMEDIATOR_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceTestUnknownFunction.hpp"
#include "SundanceTestDiscreteFunction.hpp"
#include "SundancePoint.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceMap.hpp"

namespace SundanceTesting
{
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using Sundance::Map;
  /**
   *
   */
  class TestEvalMediator : public AbstractEvalMediator
  {
  public:
    /** */
    TestEvalMediator(const Expr& fields);

    /** */
    virtual ~TestEvalMediator(){;}

    /** */
    void setEvalPoint(const Point& x) {x_=x;}

    /** */
    int numFields() const {return fields_.size();}

    /** */
    const std::string& fieldName(int i) const {return fieldNames_[i];}

    /** */
    void setFieldCoeff(int i, double A) {fields_[i].setCoeff(A);}

    /** */
    double fieldCoeff(int i) const {return fields_[i].coeff();}

    /** */
    const Map<int, int>& funcIdToFieldNumberMap() const
    {return funcIdToFieldNumberMap_;}

    /** Evaluate the given coordinate expression, putting
     * its numerical values in the given LoadableVector. */
    virtual void evalCoordExpr(const CoordExpr* expr,
                               RCP<EvalVector>& vec) const ;

    /** Evaluate the given cell diameter expression, putting
     * its numerical values in the given EvalVector. */
    virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                      RCP<EvalVector>& vec) const ;

    /** Evaluate the given cell vector expression, putting
     * its numerical values in the given EvalVector. */
    virtual void evalCellVectorExpr(const CellVectorExpr* expr,
                                      RCP<EvalVector>& vec) const ;

    /** Evaluate the given discrete function, putting
     * its numerical values in the given LoadableVector. */
    virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                         const Array<MultiIndex>& mi,
                                         Array<RCP<EvalVector> >& vec) const ;

    double evalDummyBasis(int m, const MultiIndex& mi) const ;


  private:

    Point x_;
    Map<int, int> funcIdToFieldNumberMap_;
    Array<ADField> fields_;
    Array<string> fieldNames_;
  };
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
