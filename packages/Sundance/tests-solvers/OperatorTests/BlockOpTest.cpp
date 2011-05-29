//@HEADER
// ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER


#include <cstdlib>
#include "Teuchos_GlobalMPISession.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFEpetraMatrix.hpp"
#include "TSFMatrixLaplacian1D.hpp"
#include "TSFRandomSparseMatrixBuilderDecl.hpp"
#include "TSFRandomBlockMatrixBuilderDecl.hpp"
#include "TSFCompoundTester.hpp"
#include "SundanceOut.hpp"
#include "TSFProductVectorSpaceDecl.hpp"
#include "TSFLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFProductVectorSpaceImpl.hpp"
#include "TSFRandomBlockMatrixBuilderImpl.hpp"
#endif



using namespace Teuchos;
using namespace Sundance;
using namespace TSFExtended;
using namespace TSFExtendedOps;
using Thyra::TestSpecifier;
using std::endl;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
    {
      GlobalMPISession session(&argc, &argv);
      MPIComm::world().synchronize();

      Out::os() << "go!" << std::endl;
      VectorType<double> type = new EpetraVectorType();

      Array<int> domainBlockSizes = tuple(2,3,4);
      Array<int> rangeBlockSizes = tuple(2,2);

      Array<VectorSpace<double> > domainBlocks(domainBlockSizes.size());
      Array<VectorSpace<double> > rangeBlocks(rangeBlockSizes.size());

      for (int i=0; i<domainBlocks.size(); i++)
        {
          domainBlocks[i] = type.createEvenlyPartitionedSpace(MPIComm::world(),
                                                              domainBlockSizes[i]);
        }

      for (int i=0; i<rangeBlocks.size(); i++)
        {
          rangeBlocks[i] = type.createEvenlyPartitionedSpace(MPIComm::world(),
                                                             rangeBlockSizes[i]);
        }
      
      VectorSpace<double> domain = productSpace(domainBlocks);
      VectorSpace<double> range = productSpace(rangeBlocks);

      double blockDensity = 0.75;
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;
      
      RandomBlockMatrixBuilder<double> builder(domain, range, 
        blockDensity,
        onProcDensity,
        offProcDensity,
        type);

      LinearOperator<double> A = builder.getOp();

      Out::os() << "A num block rows = " << A.numBlockRows() << std::endl;
      Out::os() << "A num block cols = " << A.numBlockCols() << std::endl;

      Vector<double> x = domain.createMember();
      Out::os() << "randomizing trial vector" << std::endl;
      Thyra::randomize(-1.0, 1.0, x.ptr().ptr());

      Array<Vector<double> > xBlock(domain.numBlocks());
      for (int i=0; i<xBlock.size(); i++)
        {
          xBlock[i] = x.getBlock(i);
        }

      Vector<double> xx = x.copy();

      

      Out::os() << "------------------------------------------------------------" << std::endl;
      Out::os() << "computing A*x..." << std::endl;
      Vector<double> y0 = A * x;
      for (int i=0; i<y0.space().numBlocks(); i++)
        {
          Out::os() << "y0[" << i << "] = " << std::endl << y0.getBlock(i) << std::endl;
        }
      

      Vector<double> y1 = range.createMember();
      Out::os() << "------------------------------------------------------------" << std::endl;
      Out::os() << "computing A*x block-by-block..." << std::endl;
      Array<Vector<double> > yBlock(range.numBlocks());
      for (int i=0; i<yBlock.size(); i++)
        {
          yBlock[i] = range.getBlock(i).createMember();
          yBlock[i].zero();
          for (int j=0; j<xBlock.size(); j++)
            {
              LinearOperator<double> Aij = A.getBlock(i,j);
              if (Aij.ptr().get() != 0)
                {
                  Out::os() << "A(" << i << ", " << j << ") = " << std::endl 
                       << Aij << std::endl;
                }
              else
                {
                  Out::os() << "A(" << i << ", " << j << ") = 0 " << std::endl;
                }
              Out::os() << "x[" << j << "] = " << std::endl << xBlock[j] << std::endl;
              if (Aij.ptr().get()==0) continue;
              yBlock[i] = yBlock[i] + Aij * xBlock[j];
            }
          y1.setBlock(i, yBlock[i]);
        }

      for (int i=0; i<y1.space().numBlocks(); i++)
        {
          Out::os() << "y1[" << i << "] = " << std::endl << y1.getBlock(i) << std::endl;
        }
      double err = (y1 - y0).norm2();
      Out::os() << "error = " << err << std::endl;

      double tol = 1.0e-13;
      if (err < tol)
        {
          Out::os() << "block op test PASSED" << std::endl;
        }
      else
        {
          stat = -1;
          Out::os() << "block op test FAILED" << std::endl;
        }
    }
  catch(std::exception& e)
    {
      stat = -1;
      Out::os() << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



