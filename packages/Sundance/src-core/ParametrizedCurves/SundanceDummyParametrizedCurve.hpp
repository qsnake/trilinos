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

#ifndef SUNDANCEDUMMYPARAMETRIZEDCURVE_H_
#define SUNDANCEDUMMYPARAMETRIZEDCURVE_H_

#include "SundanceCurveBase.hpp"

namespace Sundance
{

/** Dummy curve */
class DummyParametrizedCurve : public CurveBase
{
public:

	DummyParametrizedCurve(){;}

	virtual ~DummyParametrizedCurve(){;}

	/** */
	virtual Expr getParams() const { return Expr();}

	/** */
	virtual void getIntegrationParams(Array<double>& alphas) const {
		alphas.resize(2);
		alphas[0] = 1.0;
		alphas[1] = 1.0;
	}

	/** */
	virtual double integrationParameter(const Point& evaluationPoint) const{ return 1.0;}

	virtual double curveEquation(const Point& evaluationPoint) const { return 1.0; }

	virtual void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<Point>& result) const {;}

	virtual void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<double>& result) const {;}

	/** Return a ref count pointer to self */
	virtual RCP<CurveBase> getRcp() {return rcp(this);}

	virtual bool isCurveValid() const { return false; }


};

}
#endif /* SUNDANCEDUMMYPARAMETRIZEDCURVE_H_ */
