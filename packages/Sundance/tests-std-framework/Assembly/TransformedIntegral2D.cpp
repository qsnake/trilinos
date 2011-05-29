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


static Time& totalTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


int main(int argc, char** argv)
{
  int stat = 0;
  int verb=1;
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    int pMax = 2;
    int dim=2;

    bool isInternalBdry = false;

    Utils::setChopVal(1.0e-14);

    CellType cellType = TriangleCell;

    //       Point a = Point(1.0, 1.0);
    //       Point b = Point(1.2, 1.6);
    //       Point c = Point(0.8, 1.3);

    Point a = Point(0.0, 0.0);
    Point b = Point(1.0, 0.0);
    Point c = Point(0.0, 1.0);

    Point d = Point(0.1, 0.1);
    Point e = Point(1.0, 0.0);
    Point f = Point(0.0, 1.0);

    int nCells = 2;

    CellJacobianBatch JBatch;
    JBatch.resize(nCells, 2, 2);
    double* J = JBatch.jVals(0);
    J[0] = b[0] - a[0];
    J[1] = c[0] - a[0];
    J[2] = b[1] - a[1];
    J[3] = c[1] - a[1];

    J[4] = e[0] - d[0];
    J[5] = f[0] - d[0];
    J[6] = e[1] - d[1];
    J[7] = f[1] - d[1];


      
    Array<int> dummy;
    double coeff = 1.0;
    RCP<Array<double> > A = rcp(new Array<double>());
    RCP<Array<double> > B = rcp(new Array<double>());

    QuadratureFamily q4 = new GaussianQuadrature(4);

    int nErrors = 0;

    std::cerr << std::endl << std::endl 
         << "---------------- One-forms --------------------" 
         << std::endl << std::endl;
    for (int p=0; p<=pMax; p++)
    {
      BasisFamily P = new Lagrange(p);
      for (int dp=0; dp<=1; dp++)
      {
        if (dp > p) continue;

        int numTestDir = 1;
        if (dp==1) numTestDir = dim;
        for (int t=0; t<numTestDir; t++)
        {
          int alpha = t;
          Tabs tab;

          ParametrizedCurve curve = new DummyParametrizedCurve();
          MeshType meshType = new BasicSimplicialMeshType();
          MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
          Mesh mesh = mesher.getMesh();
          RCP<Array<int> > cellLIDs;

          RefIntegral ref(dim, cellType, dim, cellType, P, alpha, dp, q4 , isInternalBdry, curve, mesh ,verb);
          A->resize(JBatch.numCells() * ref.nNodes());
          for (int ai=0; ai<A->size(); ai++) (*A)[ai]=0.0;
          ref.transformOneForm(JBatch, JBatch, dummy, cellLIDs , coeff, A);
          std::cerr << tab << "transformed reference element" << std::endl;
          if (dp>0) std::cerr << tab << "test diff direction=" << t << std::endl;
          for (int cell=0; cell<nCells; cell++)
          {
            std::cerr << tab << "{";
            for (int r=0; r<ref.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << Utils::chop((*A)[cell*ref.nNodesTest()+r]);
            }
            std::cerr << "}" << std::endl;
          }
          QuadratureIntegral quad(dim, cellType, dim, cellType, P, alpha, dp, q4, isInternalBdry, curve, mesh, verb);
          Array<double> quadCoeff(2*quad.nQuad(), 1.0);
          B->resize(JBatch.numCells() * quad.nNodes());
          for (int ai=0; ai<B->size(); ai++) (*B)[ai]=0.0;
          quad.transformOneForm(JBatch, JBatch, dummy, cellLIDs , &(quadCoeff[0]), B);
          std::cerr << tab << "transformed quad element" << std::endl;
          if (dp>0) std::cerr << tab << "test diff direction =" << t << std::endl;
          for (int cell=0; cell<nCells; cell++)
          {
            std::cerr << tab << "{";
            for (int r=0; r<quad.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << Utils::chop((*B)[cell*ref.nNodesTest()+r]);
            }
            std::cerr << "}" << std::endl;
          }

          std::cerr << tab << "MISFIT quad-ref" << std::endl;
          std::cerr << tab << "test diff order =" << dp << std::endl;
          if (dp>0) std::cerr << tab << "test diff direction =" << t << std::endl;
          bool OK = true;
          for (int cell=0; cell<nCells; cell++)
          {
            std::cerr << tab << "{";
            for (int r=0; r<quad.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              int i = cell*ref.nNodesTest()+r;
              double err = fabs(Utils::chop((*B)[i] - (*A)[i]));
              if (err > 1.0e-14) 
              {
                OK = false;
              }
              std::cerr << err;
            }
            std::cerr << "}" << std::endl;
          }
                  
          if (!OK) 
          {
            nErrors ++;
            std::cerr << "ERROR DETECTED!!! p=" << p
                 << "  t=" << t  << std::endl;
          }
        }
      }
    }
         




    std::cerr << std::endl << std::endl 
         << "---------------- Two-forms --------------------" 
         << std::endl << std::endl;

    for (int p=0; p<=pMax; p++)
    {
      BasisFamily P = new Lagrange(p);
      for (int dp=0; dp<=1; dp++)
      {
        if (dp > p) continue;
        int numTestDir = 1;
        if (dp==1) numTestDir = dim;
        for (int q=0; q<=pMax; q++)
        {
          BasisFamily Q = new Lagrange(q);
          for (int dq=0; dq<=1; dq++)
          {
            if (dq > q) continue;
            for (int t=0; t<numTestDir; t++)
            {
              int alpha = t;
              int numUnkDir = 1;
              if (dq==1) numUnkDir = dim;
              for (int u=0; u<numUnkDir; u++)
              {
                   ParametrizedCurve curve = new DummyParametrizedCurve();
                   MeshType meshType = new BasicSimplicialMeshType();
                   MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
                   Mesh mesh = mesher.getMesh();
                   QuadratureFamily quad_1 = new GaussianQuadrature(2);
                   RCP<Array<int> > cellLIDs;

                Tabs tab;
                //                              if (p==0 || q==0 || dp==0 || dq==0 || u==1
                //  || t==1) continue;
                int beta = u;
                RefIntegral ref(dim, cellType, dim, cellType, P, alpha,
                  dp, Q, beta, dq, quad_1 , isInternalBdry, curve , mesh , verb);
                A->resize(JBatch.numCells() * ref.nNodes());
                for (int ai=0; ai<A->size(); ai++) (*A)[ai]=0.0;
                ref.transformTwoForm(JBatch, JBatch, dummy, cellLIDs , coeff, A);
                std::cerr << tab << "transformed ref element" << std::endl;
                std::cerr << tab << "test diff order = " << dp << std::endl;
                if (dp>0) std::cerr << tab << "t=dx(" << t << ")" << std::endl;
                std::cerr << tab << "unk diff order = " << dq << std::endl;
                if (dq>0) std::cerr << tab << "u=dx(" << u << ")" << std::endl;

                for (int cell=0; cell<nCells; cell++)
                {
                  std::cerr << tab << "cell=" << cell << " {";
                  for (int r=0; r<ref.nNodesTest(); r++)
                  {
                    if (r!=0) std::cerr << ", ";
                    std::cerr << "{";
                    for (int c=0; c<ref.nNodesUnk(); c++)
                    {
                      if (c!=0) std::cerr << ", ";
                      std::cerr << Utils::chop((*A)[r + ref.nNodesTest()*(c + cell*ref.nNodesUnk())]);
                    }
                    std::cerr << "}";
                  }
                  std::cerr << "}" << std::endl;
                }


                QuadratureIntegral quad(dim, cellType, dim, cellType, P, alpha,
                  dp, Q, beta, dq, q4, isInternalBdry,curve , mesh , verb);
                Array<double> quadCoeff(2*quad.nQuad(), 1.0);
                B->resize(JBatch.numCells() * quad.nNodes());
                for (int ai=0; ai<B->size(); ai++) (*B)[ai]=0.0;
                quad.transformTwoForm(JBatch, JBatch, dummy, cellLIDs , &(quadCoeff[0]), B);

                std::cerr << tab << "transformed quad element" << std::endl;
                std::cerr << tab << "test diff order = " << dp << std::endl;
                if (dp>0) std::cerr << tab << "t=dx(" << t << ")" << std::endl;
                std::cerr << tab << "unk diff order = " << dq << std::endl;
                if (dq>0) std::cerr << tab << "u=dx(" << u << ")" << std::endl;

                for (int cell=0; cell<nCells; cell++)
                {
                  std::cerr << tab << "cell=" << cell << " {";
                  for (int r=0; r<ref.nNodesTest(); r++)
                  {
                    if (r!=0) std::cerr << ", ";
                    std::cerr << "{";
                    for (int c=0; c<ref.nNodesUnk(); c++)
                    {
                      if (c!=0) std::cerr << ", ";
                      std::cerr << Utils::chop((*B)[r + ref.nNodesTest()*(c + cell*ref.nNodesUnk())]);
                    }
                    std::cerr << "}";
                  }
                  std::cerr << "}" << std::endl;   
                }

                bool OK = true;
                std::cerr << tab << "MISMATCH quad - ref" << std::endl;
                std::cerr << tab << "test diff order = " << dp << std::endl;
                if (dp>0) std::cerr << tab << "t=dx(" << t << ")" << std::endl;
                std::cerr << tab << "unk diff order = " << dq << std::endl;
                if (dq>0) std::cerr << tab << "u=dx(" << u << ")" << std::endl;

                for (int cell=0; cell<nCells; cell++)
                {
                  std::cerr << tab << "cell #" << cell << " {";
                              
                  for (int r=0; r<ref.nNodesTest(); r++)
                  {
                    if (r!=0) std::cerr << ", ";
                    std::cerr << "{";
                    for (int c=0; c<ref.nNodesUnk(); c++)
                    {
                      if (c!=0) std::cerr << ", ";
                      int i = r + ref.nNodesTest()*(c + cell*ref.nNodesUnk());
                      double err = fabs(Utils::chop((*B)[i] - (*A)[i]));
                      if (err > 1.0e-14) OK = false;
                      std::cerr << err;
                    }
                    std::cerr << "}";
                  }
                  std::cerr << "}" << std::endl;
                }
                if (!OK) 
                {
                  nErrors ++;
                  std::cerr << "ERROR DETECTED!!! p=" << p
                       << " dp=" << dp << "  t=" << t  
                       << " q=" << q << "  dq=" << dq
                       << "  u=" << u << std::endl;
                }

                std::cerr << std::endl << std::endl << std::endl << std::endl;
              }
            }
          }
        }
      }
    }

    std::cerr << "total quadrature flops: " << QuadratureIntegral::totalFlops() 
         << std::endl;
    std::cerr << "total ref integration flops: " << RefIntegral::totalFlops() 
         << std::endl;

    if (nErrors == 0)
    {
      std::cerr << "Transformed integral test PASSED" << std::endl;
    }
    else
    {
      stat = -1;
      std::cerr << "Transformed integral test FAILED" << std::endl;
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    stat = -1;
    std::cerr << "Transformed integral test FAILED" << std::endl;
    std::cerr << e.what() << std::endl;
  }

  return stat;
  
}
