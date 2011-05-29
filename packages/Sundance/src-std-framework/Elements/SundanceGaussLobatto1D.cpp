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

#include "SundanceExceptions.hpp"
#include "SundanceGaussLobatto1D.hpp"
#ifdef _MSC_VER
# include "winmath.h"
#endif

using namespace Sundance;
using namespace Teuchos;

GaussLobatto1D::GaussLobatto1D(int n) :
	nodes_(n), weights_(n)
{
	computeWeights(n, -1.0, 1.0);
}

GaussLobatto1D::GaussLobatto1D(int n, double a, double b) :
	nodes_(n), weights_(n)
{
	computeWeights(n, a, b);
}

void GaussLobatto1D::computeWeights(int n, double a, double b)
{

	TEST_FOR_EXCEPTION(n < 2, RuntimeError, "number of points=" << n
			<< " must be at least 2 for Gauss-Lobatto-Legendre quadrature!");

	int m = (n + 1) / 2;

	double xMid = (b + a) / 2.0;
	double halfWidth = (b - a) / 2.0;

	for (int i = 0; i < m; i++)
	{
		// Initial guess (Gauss-Lobatto-Chebyshev roots)
		double z = cos(M_PI * i / (n - 1));
		double p1;
		double zOld;
		double tol = 1.0e-14;
		// Find the roots of L_(n-1)' (Newton's method)
		do
		{
			p1 = 1.0;
			double p2 = 0.0;
			for (int j = 1; j < n; j++)
			{
				double p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}

			zOld = z;

			// p1' == (n-1)*(z*p1-p2)/(z*z-1) (cf. Wolfram MathWorld)
			// p1'' == ((n-1)*n*p1-2*z*p1')/(z*z-1) (Legendre differential equation,
			// neglect last summand in numerator since p1' -> 0 and abs(z)<=1)
			// This results in following loop:
			z = zOld - (z * p1 - p2) / (n * p1);
		} while (fabs(z - zOld) > tol);

		if (i == 0)
		{
			nodes_[0] = a;
			nodes_[n - 1] = b;
		}
		else
		{
			nodes_[i] = xMid - halfWidth * z;
			nodes_[n - i - 1] = xMid + halfWidth * z;
		}
		weights_[i] = 2.0 * halfWidth / ((n - 1) * n * p1 * p1);
		weights_[n - i - 1] = weights_[i];
	}
}

bool GaussLobatto1D::unitTest()
{
	std::cerr
			<< "------------------ GaussLobatto1D unit test ----------------------"
			<< std::endl;

	GaussLobatto1D q(20, 0.0, M_PI);

	double sum = 0.0;
	for (int i = 0; i < q.nPoints(); i++)
	{
		sum += q.weights()[i] * sin(q.nodes()[i]);
	}
	std::cerr << "integral of sin(x) over [0, pi] = " << sum << std::endl;
	double sumErr = fabs(sum - 2.0);
	bool sumPass = sumErr < 1.0e-10;
	std::cerr << "error = " << sumErr << std::endl;
	if (sumPass)
		std::cerr << "GaussLobatto1D sine test PASSED" << std::endl;
	else
		std::cerr << "GaussLobatto1D sine test FAILED" << std::endl;
	std::cerr
			<< "------------------ End GaussLobatto1D unit test ----------------------"
			<< std::endl;
	return sumPass;
}

