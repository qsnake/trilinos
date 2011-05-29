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

#include "SundanceBox2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Box2D::Box2D(double px, double py, double ox, double oy, double a1, double a2) :
	CurveBase(1, a1, a2), px_(px), py_(px), ox_(ox), oy_(oy)
{
}

Box2D::~Box2D()
{
}

Expr Box2D::getParams() const
{
	// return the parameters of the box
	return Expr(List(px_, py_, ox_ , oy_));
}

double Box2D::curveEquation(const Point& evalPoint) const
{
	TEST_FOR_EXCEPTION(evalPoint.dim() != 2, RuntimeError,
			"Box2D::curveEquation() evaluation point dimension must be 2");

	// calculate the distance compared to the middle point
	double distX =  fabs(px_ + 0.5*ox_ - evalPoint[0]) - 0.5*ox_;
	double distY =  fabs(py_ + 0.5*oy_ - evalPoint[1]) - 0.5*oy_;
	return (distX > distY) ? distX : distY ;
}

void Box2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	Array<double> t;
	returnIntersect(start, end, nrPoints, t);

	result.resize(nrPoints);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		result[i] = start + t[i] * (end - start);
	}
}

void Box2D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{

    // first implementation
	//TEST_FOR_EXCEPTION( true , RuntimeError, "Box2D::returnIntersect() not implemented yet");

	// calc cut in X direction
	if ( fabs(start[0] - end[0]) < 1e-6 ){
        //
	} else {
       double a = (start[1] - end[1]) / (start[0] - end[0]);
       double b = start[1] - a*start[0];
       double y1 = a*px_ + b;
       double y2 = a*(px_ + ox_) + b;
       Point p1( px_ , y1 );
       Point p2( px_ + ox_ , y2 );
       if ( (px_ >= p1[0] ) && ( px_ + ox_ <= p1[0] ) && (start[0] <= p1[0] ) && ( end[0]  >= p1[0] ) &&
    		(py_ >= p1[1] ) && ( py_ + oy_ <= p1[1] ) && (start[1] <= p1[1] ) && ( end[1]  >= p1[1] )){
    	   result.resize(1);
    	   result[0] = (p1[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
       }
       if ( (px_ >= p2[0] ) && ( px_ + ox_ <= p2[0] ) && (start[0] <= p2[0] ) && ( end[0]  >= p2[0] ) &&
    		(py_ >= p2[1] ) && ( py_ + oy_ <= p2[1] ) && (start[1] <= p2[1] ) && ( end[1]  >= p2[1] )){
    	   result.resize(2);
    	   result[1] = (p2[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
       }
	}

	// calc cut in Y direction
	if ( fabs(start[1] - end[1]) < 1e-6 ){
        //
	} else {
       double a = (start[1] - end[1]) / (start[0] - end[0]);
       double b = start[1] - a*start[0];
       double x1 = (py_ - b) * a;
       double x2 = (py_ + oy_ - b) * a;
       Point p1( x1 , py_ );
       Point p2( x2 , py_ + oy_ );
       if ( (px_ >= p1[0] ) && ( px_ + ox_ <= p1[0] ) && (start[0] <= p1[0] ) && ( end[0]  >= p1[0] ) &&
    		(py_ >= p1[1] ) && ( py_ + oy_ <= p1[1] ) && (start[1] <= p1[1] ) && ( end[1]  >= p1[1] )){
    	   result.resize(3);
    	   result[2] = (p1[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
       }
       if ( (px_ >= p2[0] ) && ( px_ + ox_ <= p2[0] ) && (start[0] <= p2[0] ) && ( end[0]  >= p2[0] ) &&
    		(py_ >= p2[1] ) && ( py_ + oy_ <= p2[1] ) && (start[1] <= p2[0] ) && ( end[1]  >= p2[1] )){
    	   result.resize(4);
    	   result[3] = (p2[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
       }
	}
}

