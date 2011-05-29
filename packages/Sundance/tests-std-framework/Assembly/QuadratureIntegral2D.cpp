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
#include "SundanceTabs.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExpr.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
//#include "SundanceFIATLagrange.hpp"
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
#include "SundanceQuadratureIntegral.hpp"
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


double chop(double x) {if (::fabs(x) < 1.0e-14) return 0.0; return x;}

int main(int argc, char** argv)
{
  
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    int pMax = 1;
    int dim=2;

    CellType cellType = TriangleCell;

    Point a = Point(0.0, 0.0);
    Point b = Point(1.0, 0.0);
    Point c = Point(0.0, 1.0);
    CellJacobianBatch JBatch;
    JBatch.resize(1, 2, 2);
    double* J = JBatch.jVals(0);
    J[0] = b[0] - a[0];
    J[1] = c[0] - a[0];
    J[2] = b[1] - a[1];
    J[3] = c[1] - a[1];


    bool isInternalBdry=false;

    /* ------ evaluate Lagrange and FIAT-Lagrange at the vertices */
    Array<Point> verts = tuple(a,b,c);
    BasisFamily lagrange = new Lagrange(1);
    BasisFamily fiatLagrange = new Lagrange(1);
      
    MultiIndex d0(0,0,0);
    MultiIndex dx(1,0,0);
    MultiIndex dy(0,1,0);

    Array<Array<Array<double> > > result;

    Array<int> dummy;

    std::cerr << "------ Evaluating bases at vertices ----------" << std::endl
         << std::endl;

    std::cerr << "Evaluating phi(vert) with FIAT-Lagrange" << std::endl;
    fiatLagrange.ptr()->refEval(cellType, verts, d0, result);
    std::cerr << "results = " << result << std::endl << std::endl;

    std::cerr << "Evaluating phi(vert) with Lagrange" << std::endl;
    lagrange.ptr()->refEval(cellType, verts, d0, result);
    std::cerr << "results = " << result << std::endl << std::endl;

    std::cerr << std::endl ;

    std::cerr << "Evaluating Dx*phi(vert) with FIAT-Lagrange" << std::endl;
    fiatLagrange.ptr()->refEval(cellType, verts, dx, result);
    std::cerr << "results = " << result << std::endl << std::endl;

    std::cerr << "Evaluating Dx*phi(vert) with Lagrange" << std::endl;
    lagrange.ptr()->refEval(cellType, verts, dx, result);
    std::cerr << "results = " << result << std::endl << std::endl;

    std::cerr << std::endl ;
      
    std::cerr << "Evaluating Dy*phi(vert) with FIAT-Lagrange" << std::endl;
    fiatLagrange.ptr()->refEval(cellType, verts, dy, result);
    std::cerr << "results = " << result << std::endl << std::endl;

    std::cerr << "Evaluating Dy*phi(vert) with Lagrange" << std::endl;
    lagrange.ptr()->refEval(cellType, verts, dy, result);
    std::cerr << "results = " << result << std::endl << std::endl;

      

    /* --------- evaluate integrals over elements ----------- */
      
    RCP<Array<double> > A = rcp(new Array<double>());
          
    QuadratureFamily quad = new GaussianQuadrature(4);
    Array<double> quadWeights;
    Array<Point> quadPts;
    quad.getPoints(cellType, quadPts, quadWeights);
    int nQuad = quadPts.size();

    Array<double> coeff(nQuad);
    for (int i=0; i<nQuad; i++) 
    {
      double s = quadPts[i][0];
      double t = quadPts[i][1];
      double x = a[0] + J[0]*s + J[1]*t;
      double y = a[1] + J[2]*s + J[3]*t;
      coeff[i] = x*y;
    }
    const double* const f = &(coeff[0]);

    std::cerr << std::endl << std::endl 
         << "---------------- One-forms --------------------" 
         << std::endl << std::endl;
    for (int p=1; p<=pMax; p++)
    {
      BasisFamily P = new Lagrange(p);
      for (int dp=0; dp<=1; dp++)
      {
        if (dp > p) continue;
        Tabs tab0;
        std::cerr << tab0 << "test function deriv order = " << dp << std::endl;
        int numTestDir = 1;
        if (dp==1) numTestDir = dim;
        for (int t=0; t<numTestDir; t++)
        {
          int alpha = t;
          Tabs tab;
          QuadratureIntegral ref(dim, cellType, dim, cellType, P, alpha, dp, quad, isInternalBdry);
          A->resize(ref.nNodesTest());
          ref.transformOneForm(JBatch, JBatch, dummy, f, A);
          std::cerr << tab << "test deriv direction =" << t << std::endl;
          std::cerr << tab << "transformed local vector: " << std::endl;
          std::cerr << tab << "{";
          for (int r=0; r<ref.nNodesTest(); r++)
          {
            if (r!=0) std::cerr << ", ";
            std::cerr << (*A)[r];
          }
          std::cerr << "}" << std::endl << std::endl;
        }
      }
    }

    std::cerr << std::endl << std::endl 
         << "---------------- Two-forms --------------------" 
         << std::endl << std::endl;
    for (int p=1; p<=pMax; p++)
    {
      BasisFamily P = new Lagrange(p);
      for (int q=1; q<=pMax; q++)
      {
        BasisFamily Q = new Lagrange(q);
        for (int dp=0; dp<=1; dp++)
        {
          if (dp > p) continue;
          Tabs tab0;
          std::cerr << tab0 << "test function deriv order = " << dp << std::endl;
          for (int dq=0; dq<=1; dq++)
          {
            if (dq > q) continue;
            Tabs tab1;
            std::cerr << tab1 
                 << "unk function deriv order = " << dq << std::endl;
            int numTestDir = 1;
            if (dp==1) numTestDir = dim;
            for (int t=0; t<numTestDir; t++)
            {
              int alpha = t;
              int numUnkDir = 1;
              if (dq==1) numUnkDir = dim;
              for (int u=0; u<numUnkDir; u++)
              {
                Tabs tab;
                int beta = u;
                QuadratureIntegral ref(dim, cellType, dim, cellType, P, alpha, 
                  dp, Q, beta, dq, quadd, isInternalBdry);
                A->resize(ref.nNodesTest()*ref.nNodesUnk());
                ref.transformTwoForm(JBatch, JBatch, dummy, f, A);

                std::cerr << tab << "test deriv direction =" << 
                  t << ", unk deriv direction =" << u << std::endl;
                std::cerr << tab << "transformed local stiffness matrix" << std::endl;
                std::cerr << tab << "{";

                for (int r=0; r<ref.nNodesTest(); r++)
                {
                  if (r!=0) std::cerr << ", ";
                  std::cerr << "{";
                  for (int c=0; c<ref.nNodesUnk(); c++)
                  {
                    if (c!=0) std::cerr << ", ";
                    std::cerr << chop((*A)[r + ref.nNodesTest()*c]);
                  }
                  std::cerr << "}";
                }
                std::cerr << "}" << std::endl << std::endl;
              }
            }
          }
        }
      }
    }
    TimeMonitor::summarize();

  }
	catch(std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}
