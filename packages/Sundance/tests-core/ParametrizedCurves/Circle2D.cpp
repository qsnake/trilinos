/*
 * Circle2D.cpp
 *
 *  Created on: Jan 12, 2010
 *      Author: benk
 */

#include "SundancePoint.hpp"
#include "SundanceCircle.hpp"
#include "SundanceParametrizedCurve.hpp"

using namespace Sundance;
using namespace Teuchos;

bool errorCheck(double v1, double v2, double tol)
{
	//std::cerr << " Doing test Value:" << v1 << " Target value:" << v2 << std::endl;
	if (fabs(v1 - v2) > tol)
	{
		std::cerr << " Test failed  Value:" << v1 << " Target value:" << v2
				<< std::endl;
		return false;
	}
	else
	{
		return true;
	}
}

int main(int argc, char** argv)
{
	// create circle
	double total_error = 0;
	bool pass = true;
	ParametrizedCurve circle = new Circle(0.5, 0.5, 0.3, 1, 0.0001);
	Point testpoint(0.0, 0.0), testpoint2(0.0, 0.0);
	int nr_points;
	Array<Point> pointVector;

	testpoint[0] = 0.5;
	testpoint[1] = 0.2;
	total_error = total_error + circle.curveEquation(testpoint);
	testpoint[0] = 0.2;
	testpoint[1] = 0.5;
	total_error = total_error + circle.curveEquation(testpoint);
	testpoint[0] = 0.2;
	testpoint[1] = 0.2;
	total_error = total_error + (circle.curveEquation(testpoint) - 0.09);
	testpoint[0] = 0.3;
	testpoint[1] = 0.3;
	pass = pass
			&& errorCheck(circle.curveEquation(testpoint) + 0.01, 0.0, 1e-4);

	testpoint[0] = 0.0;
	testpoint[1] = 0.0;
	testpoint2[0] = 1.0;
	testpoint2[1] = 1.0;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 2.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.28787, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.28787, 1e-4);
	pass = pass && errorCheck(pointVector[1][0], 0.71213, 1e-4);
	pass = pass && errorCheck(pointVector[1][1], 0.71213, 1e-4);

	testpoint[0] = 1.0;
	testpoint[1] = 1.0;
	testpoint2[0] = 0.0;
	testpoint2[1] = 0.0;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 2.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.71213, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.71213, 1e-4);
	pass = pass && errorCheck(pointVector[1][0], 0.28787, 1e-4);
	pass = pass && errorCheck(pointVector[1][1], 0.28787, 1e-4);

	testpoint[0] = 0.5;
	testpoint[1] = 0.5;
	testpoint2[0] = 0.0;
	testpoint2[1] = 0.0;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 1.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.28787, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.28787, 1e-4);

	testpoint[0] = 0.5;
	testpoint[1] = 0.5;
	testpoint2[0] = 1.0;
	testpoint2[1] = 1.0;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 1.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.71213, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.71213, 1e-4);

	testpoint[0] = 0.0;
	testpoint[1] = 0.1;
	testpoint2[0] = 1.0;
	testpoint2[1] = 1.1;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 2.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.24384, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.34384, 1e-4);
	pass = pass && errorCheck(pointVector[1][0], 0.65616, 1e-4);
	pass = pass && errorCheck(pointVector[1][1], 0.75616, 1e-4);

	testpoint[0] = 0.5;
	testpoint[1] = 0.6;
	testpoint2[0] = 1.0;
	testpoint2[1] = 1.1;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 1.0, 1e-4);
	pass = pass && errorCheck(pointVector[1][0], 0.65616, 1e-4);
	pass = pass && errorCheck(pointVector[1][1], 0.75616, 1e-4);

	testpoint[0] = 0.0;
	testpoint[1] = 0.1;
	testpoint2[0] = 0.5;
	testpoint2[1] = 0.6;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 1.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.24384, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.34384, 1e-4);

	testpoint[0] = 0.0;
	testpoint[1] = 0.1;
	testpoint2[0] = 1.0;
	testpoint2[1] = 2.1;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 2.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.20000, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.50000, 1e-4);
	pass = pass && errorCheck(pointVector[1][0], 0.32000, 1e-4);
	pass = pass && errorCheck(pointVector[1][1], 0.74000, 1e-4);

	testpoint[0] = 0.5;
	testpoint[1] = 0.0;
	testpoint2[0] = 0.5;
	testpoint2[1] = 1.0;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 2.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.5, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.2, 1e-4);
	pass = pass && errorCheck(pointVector[1][0], 0.5, 1e-4);
	pass = pass && errorCheck(pointVector[1][1], 0.8, 1e-4);

	testpoint[0] = 0.5;
	testpoint[1] = 0.4;
	testpoint2[0] = 0.5;
	testpoint2[1] = 1.0;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 1.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.5, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.8, 1e-4);

	testpoint[0] = 0.5;
	testpoint[1] = 0.1;
	testpoint2[0] = 0.5;
	testpoint2[1] = 0.7;
	circle.returnIntersectPoints(testpoint, testpoint2, nr_points, pointVector);
	pass = pass && errorCheck((double) nr_points, 1.0, 1e-4);
	pass = pass && errorCheck(pointVector[0][0], 0.5, 1e-4);
	pass = pass && errorCheck(pointVector[0][1], 0.2, 1e-4);

	// test the total error
	pass = pass && ((total_error > 1e-3) ? false : true); //error , tolerance

	if (pass)
	{
		std::cerr << "test PASSED" << std::endl;
		return 0;
	}
	else
	{
		std::cerr << "test FAILED" << std::endl;
		return -1;
	}
}
