/*
 * SundanceCurveIntagralCalc.hpp
 *
 *  Created on: Sep 2, 2010
 *      Author: benk
 */

#ifndef SUNDANCECURVEINTEGRALCALC_HPP_
#define SUNDANCECURVEINTEGRALCALC_HPP_

#include "SundanceParametrizedCurve.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceCellType.hpp"
#include "SundanceMesh.hpp"

namespace Sundance {

/** Class to compute the intersection/quadrature points of a cell with a curve in 2/3D */
class CurveIntegralCalc {
public:

	/** Empty Ctor*/
	CurveIntegralCalc();

	virtual ~CurveIntegralCalc() {;}

    /**
         [in]  maxCellType <br>
         [in]  const Array<Point>& cellPoints (the physical points of the cell) <br>
         [in]  paramCurve <br>
         [in]  Quadraturefamily , quadrature 1D for curve, or 2D for surface <br>
         [out] Array<Point>& curvePoints <br>
         [out] Array<Point>& curveDerivs  <br>
         [out] Array<Point>& curveNormals <br> */

    static void getCurveQuadPoints(CellType  maxCellType ,
    		                       int maxCellLID ,
    		                       const Mesh& mesh ,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

private :

    static void getCurveQuadPoints_line(
                                   CellType  maxCellType ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

};

}

#endif /* SUNDANCECURVEINTEGRALCALC_HPP_ */
