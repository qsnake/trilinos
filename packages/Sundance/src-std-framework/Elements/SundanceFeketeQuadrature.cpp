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

#include "SundanceFeketeBrickQuadrature.hpp"
#include "SundanceFeketeQuadQuadrature.hpp"
#include "SundanceFeketeQuadrature.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"
#include "SundanceGaussLobatto1D.hpp"
#include "SundanceTabs.hpp"
#include <stack>

using namespace Sundance;
using namespace Teuchos;

/* declare LAPACK subroutines */
extern "C"
{
/* LAPACK factorization */
void dgetrf_(const int* M, const int* N, double* A, const int* lda,
		const int* iPiv, int* info);

/* LAPACK inversion of factorized matrix */
void dgetri_(const int* n, double* a, const int* lda, const int* iPiv,
		double* work, const int* lwork, int* info);
}

FeketeQuadrature::FeketeQuadrature(int order) :
	QuadratureFamilyBase(order)
{
	_hasBasisCoeffs = false;
}

XMLObject FeketeQuadrature::toXML() const
{
	XMLObject rtn("FeketeQuadrature");
	rtn.addAttribute("order", Teuchos::toString(order()));
	return rtn;
}

void FeketeQuadrature::getLineRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	int p = order() + 3;
	p = p + (p % 2);
	int n = p / 2;

	quadPoints.resize(n);
	quadWeights.resize(n);

	GaussLobatto1D q1(n, 0.0, 1.0);

	for (int i = 0; i < n; i++)
	{
		quadWeights[i] = q1.weights()[i];
		quadPoints[i] = Point(q1.nodes()[i]);
	}
}

void FeketeQuadrature::getTriangleRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> w;

	FeketeTriangleQuadrature::getPoints(order(), w, x, y);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = 0.5 * w[i];
		quadPoints[i] = Point(x[i], y[i]);
	}
}

void FeketeQuadrature::getQuadRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> w;

	FeketeQuadQuadrature::getPoints(order(), w, x, y);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = w[i];
		quadPoints[i] = Point(x[i], y[i]);
	}
}

void FeketeQuadrature::getBrickRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> z;
	Array<double> w;

	FeketeBrickQuadrature::getPoints(order(), w, x, y, z);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = w[i];
		quadPoints[i] = Point(x[i], y[i], z[i]);
	}
}

void FeketeQuadrature::getAdaptedWeights(const CellType& cellType, int cellDim,
		int cellLID, int facetIndex, const Mesh& mesh,
		const ParametrizedCurve& globalCurve,
		Array<Point>& quadPoints, Array<double>& quadWeights,
		bool &weightsChanged) const
{
	double tol = 1.e-10;

	// Flag, if we have changed the values of the weights
	weightsChanged = false;

	// ToDo: Other cell types!
	if (cellType != TriangleCell)
		return;

	// pushForward likes cellLID in an array
	Array<int> cellLIDs(1);
	cellLIDs[0] = cellLID;

	//
	// First, let's see, if we can do it the quick way...
	//
	Array<Point> nodes;
	Array<Point> nodesref(3);
	nodesref[0] = Point(0.0, 0.0);
	nodesref[1] = Point(1.0, 0.0);
	nodesref[2] = Point(0.0, 1.0);
	mesh.pushForward(cellDim, cellLIDs, nodesref, nodes);

	// Intersection points (parameters) with the edges of triangle
	Array<double> ip01;
	Array<double> ip02;
	Array<double> ip12;
	int nr_ip01;
	int nr_ip02;
	int nr_ip12;
	globalCurve.returnIntersect(nodes[0], nodes[1], nr_ip01, ip01);
	globalCurve.returnIntersect(nodes[0], nodes[2], nr_ip02, ip02);
	globalCurve.returnIntersect(nodes[1], nodes[2], nr_ip12, ip12);

	// Cell is _not intersected_ by curve
	if (nr_ip01 + nr_ip02 + nr_ip12 == 0)
	{
		// ToDo: further intersection tests? (curve totally inside this cell...)
		double alpha = globalCurve.integrationParameter(nodes[0]);
		if (alpha != 1.0)
		{
			weightsChanged = true;
			for (int i = 0; i < quadWeights.size(); i++)
				quadWeights[i] *= alpha;
		}
		std::cout << "Nicht geschnitten!" << std::endl;
		return;
	}

	//
	// If we are here, it is an intersected cell and we have to do the big thing...
	//
	weightsChanged = true;
	std::cout << "Geschnitten!" << std::endl;
	// Number of investigated triangles
	int nTris = 1;

	// Create stack
	std::stack<Point> s;
	std::stack<double> sa;

	Array<double> originalWeights = quadWeights;
	ParametrizedCurve noCurve = new DummyParametrizedCurve();
	Array<double> integrals(originalWeights.size());
	for (int i = 0; i < integrals.size(); i++)
		integrals[i] = 0.0;

	// Put the whole triangle on stack
	s.push(Point(0.0, 0.0));
	s.push(Point(1.0, 0.0));
	s.push(Point(0.0, 1.0));
	sa.push(1.0);

	while (s.size() >= 3)
	{
		// Number of triangles investigated
		++nTris;

		std::cout << "PunkteStack " << s.size() << " AreaStack " << sa.size()
				<< std::endl;

		// Integration results
		bool refine = false;
		Array<double> results1(originalWeights.size());
		Array<double> results2(originalWeights.size());

		// Fetch the coordinates of current triangle
		double area = sa.top();
		sa.pop();
		for (int i = 2; i > -1; i--)
		{
			nodesref[i] = s.top();
			s.pop();
		}
		mesh.pushForward(cellDim, cellLIDs, nodesref, nodes);

		// Intersection points (parameters) with the edges of triangle
		globalCurve.returnIntersect(nodes[0], nodes[1], nr_ip01, ip01);
		globalCurve.returnIntersect(nodes[0], nodes[2], nr_ip02, ip02);
		globalCurve.returnIntersect(nodes[1], nodes[2], nr_ip12, ip12);

		//  Determine kind of intersection
		// 'No intersection' is handled also in first case
		Point x, xref;
		Point vec1, vec1ref;
		Point vec2, vec2ref;
		bool xInOmega = false;
		if ((nr_ip01 + nr_ip02 + nr_ip12 == 0) || (nr_ip01 == 1 && nr_ip02 == 1
				&& nr_ip12 == 0))
		{
			x = nodes[0];
			vec1 = nodes[2] - nodes[0];
			vec2 = nodes[1] - nodes[0];
			xref = nodesref[0];
			vec1ref = nodesref[2] - nodesref[0];
			vec2ref = nodesref[1] - nodesref[0];
			std::cout << "Fall 0" << std::endl;
		}
		else if (nr_ip12 == 1 && nr_ip01 == 1 && nr_ip02 == 0)
		{
			x = nodes[1];
			vec1 = nodes[0] - nodes[1];
			vec2 = nodes[2] - nodes[1];
			xref = nodesref[1];
			vec1ref = nodesref[0] - nodesref[1];
			vec2ref = nodesref[2] - nodesref[1];
			std::cout << "Fall 1" << std::endl;
		}
		else if (nr_ip02 == 1 && nr_ip12 == 1 && nr_ip01 == 0)
		{
			x = nodes[2];
			vec1 = nodes[1] - nodes[2];
			vec2 = nodes[0] - nodes[2];
			xref = nodesref[2];
			vec1ref = nodesref[1] - nodesref[2];
			vec2ref = nodesref[0] - nodesref[2];
			std::cout << "Fall 2" << std::endl;
		}
		else
		{
			std::cout << "Komische Schnitte: Verfeinern! Area: " << area
					<< std::endl;
			std::cout << "Schnitte: 01: " << nr_ip01 << " 02: " << nr_ip02
					<< " 12: " << nr_ip12 << std::endl;
			refine = true;
		}


		// If we found a case we can deal with, calculate integrals
		if (!refine)
		{
			if (globalCurve.curveEquation(x) < 0)
				xInOmega = true;
			integrateRegion(cellType, cellDim, order(), x, xref, vec1, vec2,
					vec1ref, vec2ref, globalCurve, results1);
			integrateRegion(cellType, cellDim, order() + 1, x, xref, vec1,
					vec2, vec1ref, vec2ref, globalCurve, results2);

			// Check results
			for (int i = 0; i < results2.size(); i++)
			{
				if (::fabs(results2[i] - results1[i]) > sqrt(0.5 * area * tol))
				{
					std::cout << "Zu ungenau: Verfeinern! Area: " << area
							<< std::endl;
					refine = true;
					break;
				}
			}
		}

		// If we had not found a case we can deal with or we are unsatisfied
		// with accuracy -> refinement
		if (refine)
		{
			// Divide triangle into 4 smaller ones (introduce new nodes at the
			// center of each edge
			s.push(nodesref[0]);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push((nodesref[0] + nodesref[2]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push(nodesref[1]);
			s.push((nodesref[1] + nodesref[2]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[2]) / 2);
			s.push((nodesref[1] + nodesref[2]) / 2);
			s.push(nodesref[2]);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[2]) / 2);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push((nodesref[1] + nodesref[2]) / 2);
			sa.push(area / 4);
		}

		// We found an accurate solution to that part of the triangle!
		else
		{
			if (xInOmega)
			{
				for (int i = 0; i < integrals.size(); i++)
					integrals[i] += results2[i];
				std::cout << "x in Omega!" << std::endl;
			}
			else
			{
				integrateRegion(cellType, cellDim, order() + 1, x, xref, vec1,
						vec2, vec1ref, vec2ref, noCurve, results1);
				for (int i = 0; i < integrals.size(); i++)
					integrals[i] += (results1[i] - results2[i]);
			}
		}
	} // Stack depleted here!

	if (s.size() != 0)
	{
		std::cout << "Problem mit stack size!" << std::endl;
	}
	std::cout << "Anzahl behandelter Teildreiecke: " << nTris << std::endl;

	Array<double> alphas;
	globalCurve.getIntegrationParams(alphas);
	for (int i = 0; i < quadWeights.size(); i++)
		quadWeights[i] = alphas[1] * integrals[i] + alphas[0]
				* (originalWeights[i] - integrals[i]);

	std::cout << "Gewichte: ";
	for (int i = 0; i < quadWeights.size(); i++)
		std::cout << quadWeights[i] << " ";
	std::cout << std::endl;
}

void FeketeQuadrature::integrateRegion(const CellType& cellType,
		const int cellDim, const int innerOrder, const Point& x,
		const Point& xref, const Point& vec1, const Point& vec2,
		const Point& vec1ref, const Point& vec2ref,
		const ParametrizedCurve& curve, Array<double>& integrals) const
{
	// Prepare output array
	integrals.resize((order() + 1) * (order() + 2) / 2); // ToDo: This is triangle-specific!
	for (int i = 0; i < integrals.size(); i++)
		integrals[i] = 0.0;

	// Get 'inner integration' points
	Array<double> innerQuadPts;
	Array<double> innerQuadWgts;
	int nrInnerQuadPts = innerOrder + 3;
	nrInnerQuadPts = (nrInnerQuadPts + nrInnerQuadPts % 2) / 2;
	GaussLobatto1D innerQuad(nrInnerQuadPts, 0.0, 1.0);

	innerQuadPts.resize(nrInnerQuadPts);
	innerQuadWgts.resize(nrInnerQuadPts);
	innerQuadPts = innerQuad.nodes();
	innerQuadWgts = innerQuad.weights();

	// Calculate intersection points
	for (int edgePos = 0; edgePos < nrInnerQuadPts; edgePos++)
	{
		// Direction of current integration
		Point ray = innerQuadPts[edgePos] * vec1
				+ (1.0 - innerQuadPts[edgePos]) * vec2;
		Point rayref = innerQuadPts[edgePos] * vec1ref + (1.0
				- innerQuadPts[edgePos]) * vec2ref;

		// Calculate intersection along ray
		int nrIntersect = 0;
		Array<double> t;
		Point rayEnd = x + ray;
		curve.returnIntersect(x, rayEnd, nrIntersect, t);

		// Without an intersection we integrate from x to the opposing edge
		double mu = 1.0;
		if (nrIntersect == 1)
		{
			mu = t[0];
		}
		else if (nrIntersect >= 2)
		{
			// ToDo: Verfeinern
			std::cout << "Problem: Mehr als ein Schnittpunkt mit Strahl "
					<< " fuer edgePos " << edgePos << std::endl;
			mu = t[1]; // ToDo: Falsch! Nur damit es weitergeht...
		}

		// Array for values of basis functions at inner integration points
		Array<double> innerIntegrals;
		innerIntegrals.resize(integrals.size());
		for (int i = 0; i < innerIntegrals.size(); i++)
			innerIntegrals[i] = 0.0;

		Array<Point> q;
		q.resize(nrInnerQuadPts);
		for (int rayPos = 0; rayPos < nrInnerQuadPts; rayPos++)
		{
			q[rayPos] = xref + mu * innerQuadPts[rayPos] * rayref;

			// Evaluate all basis functions at quadrature point
			Array<double> funcVals;
			evaluateAllBasisFunctions(q[rayPos], funcVals);

			// Do quadrature for all basis functions
			for (int j = 0; j < funcVals.size(); j++)
				innerIntegrals[j] += funcVals[j] * innerQuadWgts[rayPos]
						* innerQuadPts[rayPos];
		}
		for (int j = 0; j < integrals.size(); j++)
			integrals[j] += innerQuadWgts[edgePos] * mu * mu
					* innerIntegrals[j];
	}
	double det = ::fabs(vec1ref[0] * vec2ref[1] - vec2ref[0] * vec1ref[1]);
	for (int j = 0; j < integrals.size(); j++)
		integrals[j] *= det;
}

void FeketeQuadrature::evaluateAllBasisFunctions(const Point& q,
		Array<double>& result) const
{
	// Coefficients in _basisCoeffs together with PKD polynomials form
	// a Lagrange basis at Fekete points
	if (!_hasBasisCoeffs)
		computeBasisCoeffs();

	// We have nFeketePts basis functions
	int nFeketePts = sqrt(_basisCoeffs.size());
	result.resize(nFeketePts);

	// Determine values of all PKD polynomials at q
	Array<double> pkdVals;
	pkdVals.resize(nFeketePts);
	FeketeTriangleQuadrature::evalPKDpolynomials(order(), q[0], q[1],
			&(pkdVals[0]));

	// ToDo: BLAS dgemv (with matrix transposed) ?
	// For each Lagrange basis function
	for (int m = 0; m < nFeketePts; m++)
	{
		// Sum up weighted PKD polynomials
		for (int n = 0; n < nFeketePts; n++)
		{
			// Coefficients of m-th basis function are in m-th column
			result[m] += _basisCoeffs[m + n * nFeketePts] * pkdVals[n];
		}
	}
}

void FeketeQuadrature::computeBasisCoeffs() const
{
	if (_hasBasisCoeffs)
		return;

	// Get Fekete points of chosen order
	Array<Point> feketePts;

	Array<double> x;
	Array<double> y;
	Array<double> w;
	FeketeTriangleQuadrature::getPoints(order(), w, x, y);
	int nFeketePts = w.length();
	feketePts.resize(nFeketePts);
	for (int i = 0; i < nFeketePts; i++)
		feketePts[i] = Point(x[i], y[i]);

	// We construct a Lagrange basis at feketePts (there are nFeketePts of them)
	// Each basis polynomial itself is given by a linear combination of
	// nFeketePts PKD polynomials and their coefficients
	_basisCoeffs.resize(nFeketePts * nFeketePts);

	// Let's compute the coefficients:
	// Build Vandermonde matrix of PKD basis at Fekete points
	for (int n = 0; n < nFeketePts; n++)
	{
		// Set pointer to beginning of n-th row and determine values of
		// PKD polynomials at n-th Fekete point
		double* start = &(_basisCoeffs[n * nFeketePts]);
		FeketeTriangleQuadrature::evalPKDpolynomials(order(), feketePts[n][0],
				feketePts[n][1], start);
	}

	// Invert Vandermonde matrix to obtain basis coefficients:
	// LAPACK error flag, array for switched rows, work array
	int lapack_err = 0;
	Array<int> pivot;
	pivot.resize(nFeketePts);
	Array<double> work;
	work.resize(1);
	int lwork = -1;

	// LU factorization
	::dgetrf_(&nFeketePts, &nFeketePts, &(_basisCoeffs[0]), &nFeketePts,
			&(pivot[0]), &lapack_err);

	TEST_FOR_EXCEPTION(
			lapack_err != 0,
			RuntimeError,
			"FeketeQuadrature::computeBasisCoeffs(): factorization of generalized Vandermonde failed");

	// Determine work array size
	::dgetri_(&nFeketePts, &(_basisCoeffs[0]), &nFeketePts, &(pivot[0]),
			&(work[0]), &lwork, &lapack_err);
	lwork = (int) work[0];
	work.resize(lwork);

	// Inversion of factorized matrix
	::dgetri_(&nFeketePts, &(_basisCoeffs[0]), &nFeketePts, &(pivot[0]),
			&(work[0]), &lwork, &lapack_err);

	TEST_FOR_EXCEPTION(
			lapack_err != 0,
			RuntimeError,
			"FeketeQuadrature::computeBasisCoeffs(): inversion of generalized Vandermonde failed");

	_hasBasisCoeffs = true;
}
