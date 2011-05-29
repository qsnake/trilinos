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

#ifndef SUNDANCE_CURVEEVALMEDIATOR_H
#define SUNDANCE_CURVEEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * 
 */
class CurveEvalMediator : public StdFwkEvalMediator
{
public:


  /** */
  CurveEvalMediator(const Mesh& mesh,
	const ParametrizedCurve& paramcurve ,
    int cellDim,
    const QuadratureFamily& quad);


  /** */
  virtual ~CurveEvalMediator(){;}

      
  /** Evaluate the given coordinate expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCoordExpr(const CoordExpr* expr,
    RCP<EvalVector>& vec) const ;
      
  /** Evaluate the given discrete function, putting
   * its numerical values in the given EvalVector. */
  virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
    const Array<MultiIndex>& mi,
    Array<RCP<EvalVector> >& vec) const ;

  /** Evaluate the given cell diameter expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
    RCP<EvalVector>& vec) const ;

  /** Evaluates one component of the normal vector to a given parameterized curve
   * i.e. x,y or z component of that vector in 3D <br>
   * , this method is only in the CurveEvalMediator class implemented */
  virtual void evalCurveNormExpr(const CurveNormExpr* expr,
    RCP<EvalVector>& vec) const ;

  /** Evaluate the given cell vector expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellVectorExpr(const CellVectorExpr* expr,
    RCP<EvalVector>& vec) const ;
            
  /** */
  virtual void setCellType(const CellType& cellType,
    const CellType& maxCellType,
    bool isInternalBdry) ;

  /** */
  virtual void print(std::ostream& os) const ;

  /** */
  RCP<Array<Array<Array<double> > > > getFacetRefBasisVals(const BasisFamily& basis, int diffOrder) const ;

  /** */
  int numQuadPts(const CellType& cellType) const ;

  static double& totalFlops() {static double rtn = 0; return rtn;}

      

  static void addFlops(const double& flops) {totalFlops() += flops;}

  /**
   * Return the number of different cases for which reference
   * basis functions must be evaluated. If we're on maximal cells,
   * this will be one. If we're on lower-dimensional cells, this will
   * be equal to the number of cellDim-dimensional facets of the maximal
   * cells.
   */
  int numEvaluationCases() const {return numEvaluationCases_;}


  /** */
  static Time& coordEvaluationTimer() ;

private:

  /** */
  int numEvaluationCases_;

  /** */
  QuadratureFamily quad_;

  /** */
  int numQuadPtsForMaxCell_;

  /** */
  CellType maxCellType_;

  /** */
  CellType curveCellType_;

  /** */
  ParametrizedCurve paramcurve_;
      
};
}


#endif /* SUNDANCE_CURVEEVALMEDIATOR_H */
