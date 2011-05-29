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

#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureFamily.hpp"


using namespace Teuchos;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;


int iPow(int base, int exponent)
{
  int rtn = 1;
  for (int i=0; i<exponent; i++)
    {
      rtn *= base;
    }
  return rtn;
}

namespace Sundance {
void checkbasis( BasisFamily &b1 , BasisFamily &b2 )
{
  int maxDim=3;
  double tol = 1.0e-13;
  int maxDiffOrder = 0;
  int numErrors = 0;
  QuadratureFamily quad = new GaussianQuadrature(4);
  
  for (int spatialDim=1; spatialDim<=maxDim; spatialDim++) {
    std::cerr << "\t" << "spatial dimension =" << spatialDim << std::endl;
    for (int cellDim=0; cellDim<=spatialDim; cellDim++) { 
      std::cerr << "\t\t" << "cell dimension =" << cellDim << std::endl;
      CellType cellType;
      if (cellDim==0) cellType=PointCell;
      if (cellDim==1) cellType=LineCell;
      if (cellDim==2) cellType=TriangleCell;
      if (cellDim==3) cellType=TetCell;
      
      Array<Point> qPts;
      Array<double> qWts;
      quad.getPoints(cellType, qPts, qWts);
      
      for (int d=0; d<=maxDiffOrder; d++) {
	if (cellDim==0 && d>0) continue;
	cerr << "\t\t\t" << "differentiation order = " << d << std::endl;
	for (int dir=0; dir<iPow(cellDim, d); dir++) {
	  std::cerr << "\t\t\t\t" << "direction = " << dir << std::endl;
	  MultiIndex mi;
	  mi[dir]=d;
	  Array<Array<double> > values1;
	  Array<Array<double> > values2;
	  std::cerr << "\t\t\t\t" << "computing basis1...";
	  b1.ptr()->refEval(spatialDim, cellType, qPts, mi, values1);
	  std::cerr << "done" << std::endl << "\t\t\t\t" << "computing basis2...";
	  b2.ptr()->refEval(spatialDim, cellType, qPts, mi, values2);
	  std::cerr << "done" << std::endl;
	  int nNodes1 = b1.ptr()->nNodes(spatialDim, cellType);
	  int nNodes2 = b2.ptr()->nNodes(spatialDim, cellType);
	  std::cerr << "\t\t\t\t" << "num nodes: basis1=" << nNodes1
	       << " basis2=" << nNodes2 << std::endl;
	  if (nNodes1 != nNodes2) {	
	    std::cerr << "******** ERROR: node counts should be equal" << std::endl;
	    numErrors++;
	    continue;
	  }
	  if (values1.size() != values2.size()) {
	    std::cerr << "******** ERROR: value array outer sizes should be equal" << std::endl;
	    numErrors++;
	    continue;
	  }
	  if (values1.size() != qPts.size()) {
	    std::cerr << "******** ERROR: value array outer size should be equal to number of quad points" << std::endl;
	    numErrors++;
	    continue;
	  }
	  for (int q=0; q<qPts.length(); q++) {
	    if (values1[q].length() != nNodes1) {
	      std::cerr << "******** ERROR: value array inner size should be equal to number of nodes" << std::endl;
	      numErrors++;
	      continue;
	    }
	    std::cerr << "\t\t\t\t\t" << "quad point q=" << q << " pt=" << qPts[q]
		 << std::endl;
	    for (int n=0; n<nNodes1; n++) {
	      std::cerr << "\t\t\t\t\t\t" << "node n=" << n << " phi1="
		   << values1[q][n]
		   << " phi2=" << values2[q][n]
		   << " |phi1-phi2|=" << fabs(values1[q][n]-values2[q][n])
		   << std::endl;
	      if (fabs(values1[q][n]-values2[q][n]) > tol) {
		cout << "ERROR" << std::endl; numErrors++;
	      }
	    }
	  }
	}
      }
    }
  }    
  std::cerr << std::endl << std::endl << "Summary: detected " << numErrors << " errors " << std::endl;
}
}
