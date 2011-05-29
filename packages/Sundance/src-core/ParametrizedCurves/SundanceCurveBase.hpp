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


#ifndef SUNDANCECURVEBASE_H_
#define SUNDANCECURVEBASE_H_

#include "SundanceExpr.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"
#include "SundanceHandleable.hpp"

namespace Sundance
{

/**
 * The base class for parameterized curves<br>
 * This class is mean to used in the Adaptive Cell Integration. <br>
 * This class optionally can passed to one integration and the curves defines
 * where the integral is defined, and where it is not defined. <br>
 * The class defines the basic interface methods, which are needed in the
 * Adaptive Cell Integration.
 */

class CurveBase : public Sundance::Handleable<CurveBase>
{
public:

	CurveBase(){;}

	/** Base constructor
	 * @param dim , the dimension of the curve, where it is defined
	 * @param alpha1 , the integration coefficient for inside the domain
	 * @param alpha2 , the integration coefficient for outside the domain*/
	CurveBase(int dim , double alpha1 , double alpha2):_dimension(dim),
	_alpha1(alpha1),_alpha2(alpha2) {
	}

	virtual ~CurveBase(){;}

	/** This function should be used later in the optimization process
	 * TODO: This interface should be extended later
	 * @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const = 0 ;

	/**
	 * Return the integration parameters
	 * @param alphas, integration parameters, first one corresponds to Equation > 0
	 */
	void getIntegrationParams(Array<double>& alphas) const {
		alphas.resize(2);
		alphas[0] = _alpha1;
		alphas[1] = _alpha2;
	}

	/** return the integration coefficient for the inner */
	double getAlpha1() const { return _alpha1; }

	/** return the integration coefficient for the outer */
	double getAlpha2() const { return _alpha2; }

	/** return the dimension of the curve , in which it is defined */
	int getCurveDim() const { return _dimension; }

	/**
	 * @param evaluationPoint the point where we want the alpha integration parameter <br>
	 * The parameter alpha is used for the Adaptive Cell Integration
	 * @return double the alpha parameter for the integration */
	double integrationParameter(const Point& evaluationPoint) const{// here we return the integration parameter
		// this function can be overwritten
		if (curveEquation(evaluationPoint) > 0 )
			return _alpha1;
		else
			return _alpha2;
	}

	/**
	 * This function should be implemented
	 * @param evaluationPoint the point where we want the alpha integration parameter <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation(const Point& evaluationPoint) const = 0 ;

	/**
	 * This function is important for nonstructural mesh integration.<br>
	 * The function returns the intersection points with a given line (in 2D and in 3D)<br>
	 * The line is defined with a starting point and an end point
	 * @param startEdgePoint , the start point of the line
	 * @param endEdgePoint , the end point of the line
	 * @param nrPoints , number of resulted (intersected) points
	 * @param result , the resulted points of intersections */
	virtual void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<Point>& result) const = 0 ;

	/**
	 * As above, but, instead of coordinates, returns intersection values t in [0,1]
	 * and the line is defined by "start + t * (end-start)"
	 */
	virtual void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
				                      int& nrPoints ,Array<double>& result) const = 0 ;

	/** In some constructors we need to pass a Parametrized curve even if there is any
	 * so we need sometimes to pass a dummy curve, and this function shows which curve is dummy and which is not */
	virtual bool isCurveValid() const = 0;

	/** function which shows if some integral schould be calculated along the curve <br>
	 * default this returns false , only one wrapper class should overwrite this and return false */
	virtual bool isCurveIntegral() const { return false; }

protected:

	/** The dimension where the curve is defined */
	int _dimension;

	/** The integration parameter when the curve equation is positive*/
	double _alpha1;

	/** The integration parameter when the curve equation is negative*/
	double _alpha2;
};

} // end from namespace
#endif /* SUNDANCECURVEBASE_H_ */
