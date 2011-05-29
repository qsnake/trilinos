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

#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceRefIntegral.hpp"
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

int iPow(int base, int exponent)
{
  int rtn = 1;
  for (int i=0; i<exponent; i++)
  {
    rtn *= base;
  }
  return rtn;
}

int main(int argc, char** argv)
{
  int stat = 0;

  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());
    Tabs tab0;

    int pMax = 3;
    int maxDim=2;
    double tol = 1.0e-13;
    int maxDiffOrder = 0;
    int numErrors = 0;
      
    QuadratureFamily quad = new GaussianQuadrature(4);

    for (int p=1; p<=pMax; p++)
    {
      Tabs tab1;
      std::cerr << tab1 << "Polynomial order p=" << p << std::endl;
      for (int spatialDim=0; spatialDim<=maxDim; spatialDim++)
      {
        Tabs tab2;
        std::cerr << tab2 << "spatial dimension =" << spatialDim << std::endl;
        for (int cellDim=0; cellDim<=spatialDim; cellDim++)
        {   
          Tabs tab3;
          std::cerr << tab3 << "cell dimension =" << cellDim << std::endl;
          CellType cellType;
          if (cellDim==0) cellType=PointCell;
          if (cellDim==1) cellType=LineCell;
          if (cellDim==2) cellType=TriangleCell;
          if (cellDim==3) cellType=TetCell;
                  
          Array<Point> qPts;
          Array<double> qWts;
          quad.getPoints(cellType, qPts, qWts);

          BasisFamily b1 = new Lagrange(p);
          BasisFamily b2 = new Lagrange(p);

          for (int d=0; d<=maxDiffOrder; d++)
          {
            if (cellDim==0 && d>0) continue;
            Tabs tab4;
            std::cerr << tab4 << "differentiation order = " << d << std::endl;
                      
            for (int dir=0; dir<iPow(cellDim, d); dir++)
            {
              Tabs tab5;
              std::cerr << tab5 << "direction = " << dir << std::endl;
              MultiIndex mi;
              mi[dir]=d;
              Array<Array<Array<double> > > values1;
              Array<Array<Array<double> > > values2;
              std::cerr << tab5 << "computing basis1...";
              b1.ptr()->refEval(cellType, qPts, mi, values1);
              std::cerr << "done" << std::endl;
              std::cerr << tab5 << "computing basis2...";
              b2.ptr()->refEval(cellType, qPts, mi, values2);
              std::cerr << "done" << std::endl;
              int nNodes1 = b1.ptr()->nReferenceDOFsWithFacets(cellType, cellType);
              int nNodes2 = b2.ptr()->nReferenceDOFsWithFacets(cellType, cellType);
              std::cerr << tab5 << "num nodes: basis1=" << nNodes1
                   << " basis2=" << nNodes2 << std::endl;
              if (nNodes1 != nNodes2) 
              {
                std::cerr << "******** ERROR: node counts should be equal" << std::endl;
                numErrors++;
                continue;
              }
              if (values1.size() != values2.size())
              {
                std::cerr << "******** ERROR: value array outer sizes should be equal" << std::endl;
                numErrors++;
                continue;
              }
              if (values1[0].size() != qPts.size())
              {
                std::cerr << "******** ERROR: value array outer size should be equal to number of quad points" << std::endl;
                numErrors++;
                continue;
              }
              for (int q=0; q<qPts.length(); q++)
              {
                if (values1[0][q].length() != nNodes1)
                {
                  std::cerr << "******** ERROR: value array inner size should be equal to number of nodes" << std::endl;
                  numErrors++;
                  continue;
                }
                Tabs tab6;
                std::cerr << tab6 << "quad point q=" << q << " pt=" << qPts[q]
                     << std::endl;
                for (int n=0; n<nNodes1; n++)
                {
                  Tabs tab7;
                  std::cerr << tab7 << "node n=" << n << " phi1=" 
                       << values1[0][q][n] 
                       << " phi2=" << values2[0][q][n] 
                       << " |phi1-phi2|=" << fabs(values1[0][q][n]-values2[0][q][n]) 
                       << std::endl;
                  if (fabs(values1[0][q][n]-values2[0][q][n]) > tol) { cout << "ERROR" << std::endl; numErrors++; }
                }
              }
            }
          }
        }
      }
    }

    std::cerr << std::endl << std::endl << "Summary: detected " << numErrors << " errors " << std::endl;
      
    if (numErrors == 0)
    {
      std::cerr << "BasisCheck PASSED" << std::endl;
    }
    else
    {
      std::cerr << "BasisCheck FAILED" << std::endl;
      stat = -1;
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    std::cerr << "BasisCheck FAILED" << std::endl;
    stat = -1;
  }

  return stat;
}
