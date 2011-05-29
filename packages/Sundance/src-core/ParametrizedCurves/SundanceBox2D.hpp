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

#ifndef SUNDANCEBOX2D_H_
#define SUNDANCEBOX2D_H_

#include "SundanceCurveBase.hpp"

namespace Sundance
{

/** This class models a rectangle (Box) in 2D.<br>
 * This is one of the simplest curves */
class Box2D: public Sundance::CurveBase
{
public:

	/** The constructor of a 2D
	 * @param px the lowest X coordinate of the box
	 * @param py the lowest Y coordinate of the box
	 * @param ox offset of the box in X direction
	 * @param oy offset of the box in Y direction
	 * @param a1 alpha1 coefficient for FCM
	 * @param a2 alpha2 coefficient for FCM */
	Box2D(double px, double py, double ox, double oy, double a1, double a2);

	virtual ~Box2D();

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

	/**
	 * This function should be implemented
	 * @param evalPoint point where the box equation is evaluated <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation(const Point& evalPoint) const;

	/**
	 * This function is important for nonstructural mesh integration.<br>
	 * The function returns the intersection points with a given line (only 2D at the moment)<br>
	 * The line is defined with a starting point and an end point<br>
	 * The resulting points should be between the two points
	 * @param start, the start point of the line
	 * @param end, the end point of the line
	 * @param nrPoints , number of resulted (intersected) points
	 * @param result , the resulted points of intersections */
	virtual void returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
			Array<Point>& result) const;

	/**
	 * As above, but, instead of coordinates, returns intersection values t in [0,1]
	 * and the line is defined by "start + t * (end-start)"
	 */
	virtual void returnIntersect(const Point& start, const Point& end, int& nrPoints,
			Array<double>& result) const;

	/** Return a ref count pointer to self */
	virtual RCP<CurveBase> getRcp()
	{
		return rcp(this);
	}

	/** This curve is a real curve */
	virtual bool isCurveValid() const
	{
		return true;
	}

private:

	/** */
	double px_;

	/** */
	double py_;

	/** */
	double ox_;

	/** */
	double oy_;

};

}
#endif /* SUNDANCEBOX2D_H_ */
