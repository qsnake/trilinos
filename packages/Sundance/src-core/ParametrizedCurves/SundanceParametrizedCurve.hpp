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

#ifndef SUNDANCEPARAMETRIZEDCURVE_H_
#define SUNDANCEPARAMETRIZEDCURVE_H_

#include "SundanceHandle.hpp"
#include "SundanceDummyParametrizedCurve.hpp"
#include "SundanceCurveBase.hpp"

/** This specifies the different filter modes
 *  If we have a closed curve (see Circle) than we inside that curve the "curveEquation"
 *  function should be negative <br>
 *  On the outside the curveEquation function should return a positive value. */
enum CurveCellFilterMode {Outside_Curve = 1 ,
	                      Inside_Curve = -1,
	                      On_Curve = 0 };

namespace Sundance
{

/**
 * This class a "handle" to the parameterized curves, which might be used
 * in Adaptive Cell Integration.<br> e.g : the curve specifies
 * which is my computational domain and which is not.<br>
 * But this might be used for other purposes as well.
 *
 */
class ParametrizedCurve : public Sundance::Handle<CurveBase> {
public:

	/* */
	/*HANDLE_CTORS(ParametrizedCurve, CurveBase);*/

    /** Construct an empty ParametrizedCurve type object */
	ParametrizedCurve();

    /** Construct from a raw pointer to a ParametrizedCurve type subtype */
	ParametrizedCurve(Sundance::Handleable<CurveBase>* rawPtr);

    /** Construct from a smart pointer to a ParametrizedCurve type subtype */
	ParametrizedCurve(const RCP<CurveBase>& smartPtr);

	virtual ~ParametrizedCurve();

	/** Returns the parameters of the curve*/
	Expr getParams() const {
		return ptr()->getParams();
	}

	/** List integration parameters for the FCM method*/
	void getIntegrationParams(Array<double>& alphas) const {
		return ptr()->getIntegrationParams(alphas);
	}

	/** return the integration coefficient for the inner */
	double getAlpha1() const { return ptr()->getAlpha1(); }

	/** return the integration coefficient for the outer */
	double getAlpha2() const { return ptr()->getAlpha2(); }

	/** return the dimension of the curve , in which it is defined */
	int getCurveDim() const { return ptr()->getCurveDim(); }

	/** Returns the integration parameter for the FCM method*/
	double integrationParameter(const Point& evaluationPoint) const {
		return ptr()->integrationParameter(evaluationPoint);
	}

	/** The curve(2D, 3D curve, 3D surface) equation */
	double curveEquation(const Point& evaluationPoint) const {
		return ptr()->curveEquation(evaluationPoint);
	}

	/** Returns the points of intersection of the curve with the line defined by the points <br>
	 *  Only intersection between start and end point will be stated */
	void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
	                      int& nrPoints ,Array<Point>& result) const {
		ptr()->returnIntersectPoints( startEdgePoint , endEdgePoint , nrPoints , result);
	}

	/**
	 *  Same as above, but intersection parameters t will be returned where
	 *  the line is defined by start+t*(end-start). Consequently all t's will be
	 *  in the range [0,1]
	 */
	void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
		                      int& nrPoints ,Array<double>& result) const {
			ptr()->returnIntersect( startEdgePoint , endEdgePoint , nrPoints , result);
		}

	/** Shows if the curve is a valid curve*/
	inline bool isCurveValid() const { return ptr()->isCurveValid(); }

	/** function which shows if some integral schould be calculated along the curve, <br>
	 * True means we have to calculate a curve/surf integral */
	inline bool isCurveIntegral() const { return ptr()->isCurveIntegral(); }

	/** Static class returns the dummy class when this is needed */
	static const ParametrizedCurve& returnDummyCurve() {return handle_; }

	/** The operator to determine the uniqueness of one curve
	 *  compare to another ParametrizedCurve, used for placement in STL containers*/
	bool operator<(const ParametrizedCurve& other) const { return (myID_ < other.myID_)?true:false;}

	/** return the ID of the curve */
	int myID() const { return myID_; }

private:

	/** The ID to identify the ParamCurve uniquely*/
	int myID_;

	/* The dummy curve which does not do anything */
	static const ParametrizedCurve handle_;

	/** used to give global IDs*/
	static int globalId_;
};

}
#endif /* SUNDANCEPARAMETRIZEDCURVE_H_ */
