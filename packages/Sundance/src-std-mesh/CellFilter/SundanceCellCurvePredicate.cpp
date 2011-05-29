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

/*
 * SundanceCellCurvePredicate.cpp
 *
 *  Created on: Feb 19, 2010
 *      Author: benk
 */

#include "SundanceCellCurvePredicate.hpp"
#include "SundanceStdMathOps.hpp"

using namespace Sundance;

#define SIGN(X) ((X>0.0)?1:-1)

bool CellCurvePredicate::lessThan(const CellPredicateBase* other) const
{
  const CellCurvePredicate* S = dynamic_cast<const CellCurvePredicate*>(other);

  TEST_FOR_EXCEPTION( S== 0,
                     InternalError,
                     "argument " << other->toXML()
                     << " to CellCurvePredicate::lessThan() should be "
                     "a CellCurvePredicate pointer.");
  // This comparison is IMPORTANT, to determine RQCs
  return OrderedPair<ParametrizedCurve, int>(curve_, (int)filterMode_)
      < OrderedPair<ParametrizedCurve, int>(S->curve_, (int)S->filterMode_);
}

void CellCurvePredicate::testBatch(const Array<int>& cellLID,
                                        Array<int>& results) const
{
  results.resize(cellLID.size());

  if (cellDim()==0)
    {
	  switch (filterMode_){
	  case Outside_Curve:
	      for (int i=0; i<cellLID.size(); i++)
	    	  if ( curve_.curveEquation( mesh().nodePosition(cellLID[i])) > 0.0 )
	              results[i] = true;
	    	  else
	    		  results[i] = false;
		  break;
	  case Inside_Curve:
	      for (int i=0; i<cellLID.size(); i++)
	    	  if ( curve_.curveEquation( mesh().nodePosition(cellLID[i])) < 0.0 )
	              results[i] = true;
	    	  else
	    		  results[i] = false;
		  break;
	  case On_Curve:
	      for (int i=0; i<cellLID.size(); i++)
	    	  if ( fabs(curve_.curveEquation( mesh().nodePosition(cellLID[i]))) < 1e-8 )
	              results[i] = true;
	    	  else
	    		  results[i] = false;
		  break;
	  }
    }
  else
    {
      Array<int> facetLIDs;
      Array<int> facetSigns;
      int nf = mesh().numFacets(cellDim(), cellLID[0], 0);
      mesh().getFacetLIDs(cellDim(), cellLID, 0, facetLIDs, facetSigns);
      for (int c=0; c<cellLID.size(); c++)
        {
          results[c] = true;
          if (filterMode_ == On_Curve) results[c] = false;
          int curve_sign = 0;
          for (int f=0; f<nf; f++)
            {
        	  int fLID = facetLIDs[c*nf + f];
        	  switch (filterMode_){
        	  case Outside_Curve:
        		  if ( curve_.curveEquation( mesh().nodePosition(fLID) ) <= 0.0 ) {
        		     results[c] = false;
        		     continue;
        		  }
        		  break;
        	  case Inside_Curve:
        		  if ( curve_.curveEquation( mesh().nodePosition(fLID) ) >= 0.0 ) {
        		     results[c] = false;
        		     continue;
        		  }
        		  break;
        	  case On_Curve:
        		  if (f == 0){
        			  curve_sign = SIGN(curve_.curveEquation( mesh().nodePosition(fLID)));
        		  } else {
        			  if ( curve_sign != SIGN(curve_.curveEquation( mesh().nodePosition(fLID))) ){
             		     results[c] = true;
             		     continue;
        			  }
        		  }
        		  break;
        	  }
            } // from for loop
        }
    }
}

XMLObject CellCurvePredicate::toXML() const
{
  XMLObject rtn("SundanceCellCurvePredicate");
  return rtn;
}
