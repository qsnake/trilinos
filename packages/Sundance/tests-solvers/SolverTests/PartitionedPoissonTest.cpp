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

#include "Teuchos_GlobalMPISession.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFLinearCombinationDecl.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFInverseOperatorDecl.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "TSFAztecSolver.hpp"
#include "TSFMatrixLaplacian1D.hpp"
#include "TSFLinearSolverBuilder.hpp"
#include "TSFBCPartitionedVSBuilder.hpp"
#include "TSFPartitionedMatrixFactory.hpp"
#include "TSFLoadableMatrix.hpp"
#include "SundancePathUtils.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "TSFPartitionedToMonolithicConverter.hpp"
#include <set>

#include "TSFLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"

#endif

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;
using std::set;

int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
  {
    GlobalMPISession session(&argc, &argv);

    int np = MPIComm::world().getNProc();

    if (np > 1)
    {
      cout << "parallel partioned poisson test INACTIVE" << std::endl;
    }
    else
    {

#define INACTIVE
#ifdef INACTIVE
      cout << "test INACTIVE" << std::endl;
#else
      int debugWait = 0;
      if (debugWait)
      {
        int wait=1;
        int pid = getpid();
        std::cerr << "PID=" << pid << std::endl;
        std::string myCommandName=argv[0];
        std::string debugCmd = "ddd --gdb -x ~/.gdbinit " + myCommandName 
          + " " + Teuchos::toString(pid) + " &";
        std::cerr << "launching " << debugCmd << std::endl;
        system(debugCmd.c_str());
        while (wait) {;}
      }

      Teuchos::RCP<Teuchos::FancyOStream>
        out = Teuchos::VerboseObjectBase::getDefaultOStream();


      const Teuchos::EVerbosityLevel
        verbLevel = Teuchos::VERB_EXTREME;

      int rank = MPIComm::world().getRank();
      int nProc = MPIComm::world().getNProc();

      MPIComm::world().synchronize();

      VectorType<double> type = new EpetraVectorType();

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(Sundance::searchForFile("poissonParams.xml"));
#else
      ParameterXMLFileReader reader("poissonParams.xml");
#endif

      ParameterList solverParams = reader.getParameters();

      /* create the vector spaces */
      int nx = solverParams.get<int>("nLocal");
      int nLocalRows = 3 * nx;
      int nTotalRows = nProc * nLocalRows;
      int lowestLocalRow = nLocalRows * rank;
      RCP<Array<int> > isBCRow = rcp(new Array<int>(nLocalRows, 0));
      RCP<Array<int> > isBCCol = rcp(new Array<int>(nLocalRows, 0));
      RCP<set<int> > remoteBCCols = rcp(new set<int>());

      for (int i=0; i<nx; i++)
      {
        (*isBCRow)[3*i] = 1;
        (*isBCRow)[3*i+2] = 1;
        (*isBCCol)[3*i] = 1;
        (*isBCCol)[3*i+2] = 1;
        if (i==0) 
        {
          if (rank==0)
          {
            (*isBCRow)[3*i+1] = 1;
            (*isBCCol)[3*i+1] = 1;
          }
          else
          {
            remoteBCCols->insert(lowestLocalRow + 3*(i-1));
            remoteBCCols->insert(lowestLocalRow + 3*(i-1)+2);
          }
        } 
        if (i==nx-1) 
        {
          if (rank==nProc-1)
          {
            (*isBCRow)[3*i+1] = 1;
            (*isBCCol)[3*i+1] = 1;
          }
          else
          {
            remoteBCCols->insert(lowestLocalRow + 3*(i+1));
            remoteBCCols->insert(lowestLocalRow + 3*(i+1)+2);
          }
        } 
      }

      MPIComm::world().synchronize();
      VectorSpace<double> range = buildPartitionedSpace(
        nTotalRows,
        lowestLocalRow,
        nLocalRows,
        *isBCRow,
        type,
        type,
        MPIComm::world());

      VectorSpace<double> domain = buildPartitionedSpace(
        nTotalRows,
        lowestLocalRow,
        nLocalRows,
        *isBCCol,
        type,
        type,
        MPIComm::world());

      RCP<MatrixFactory<double> > mf 
        = rcp(new PartitionedMatrixFactory(domain, lowestLocalRow,
            isBCCol, remoteBCCols, type,
            range, lowestLocalRow, isBCRow, type));

      IncrementallyConfigurableMatrixFactory* icmf 
        = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(mf.get());
      for (int i=0; i<nx; i++)
      {
        for (int j=0; j<3; j++)
        {
          int row = lowestLocalRow + 3*i+j;
          Array<int> colIndices;
          if ((j==0 || j==2) 
            || (rank==0 && i==0) || (rank==nProc-1 && i==nx-1))
          {
            colIndices = tuple(row);
          }
          else
          {
            colIndices = tuple(row-4, row-3, row-2, row, row+2, row+3, row+4);
          }
          icmf->initializeNonzerosInRow(row, colIndices.size(),
            &(colIndices[0]));
        }
      }
      icmf->finalize();

      LinearOperator<double> A = mf->createMatrix();

      LoadableMatrix<double>* mat = dynamic_cast<LoadableMatrix<double>*>(A.ptr().get());

      /* fill in with the Laplacian operator */
      for (int i=0; i<nx; i++)
      {
        for (int j=0; j<3; j++)
        {
          int row = lowestLocalRow + 3*i+j;
          Array<int> colIndices;
          Array<double> colVals;
          if ((j==0 || j==2) 
            || (rank==0 && i==0) || (rank==nProc-1 && i==nx-1))
          {
            colIndices = tuple(row);
            colVals = tuple(1.0);
          }
          else
          {
            colIndices = tuple(row-4, row-3, row-2, row, row+2, row+3, row+4);
            colVals = tuple(1.0, -1.0, 1.0, 2.0, 1.0, -1.0, 1.0);
          }
          mat->addToRow(row, colIndices.size(), 
            &(colIndices[0]), &(colVals[0]));
        }
      }

#ifdef LOUD
      cout << "A = " << A << std::endl;
      for (int br=0; br<A.numBlockRows(); br++)
      {
        for (int bc=0; bc<A.numBlockCols(); bc++)
        {
          cout << "A[" << br << ", " << bc << "]=" << std::endl << A.getBlock(br,bc) << std::endl;
        }
      }
#endif


      Vector<double> x = A.domain().createMember();
      Thyra::randomize(-ST::one(),+ST::one(),x.ptr().ptr());

      Vector<double> b = A*x;

      Vector<double> ans = A.range().createMember();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      LinearOperator<double> A00 = A.getBlock(0,0);
      LinearOperator<double> A01 = A.getBlock(0,1);
      LinearOperator<double> A11 = A.getBlock(1,1);

      LinearOperator<double> A11Inv = A11.inverse(solver);
      LinearOperator<double> A00Inv = A00.inverse(solver);

      Vector<double> b0 = b.getBlock(0);
      Vector<double> b1 = b.getBlock(1);

      Vector<double> a1 = A11Inv*b1;

      Vector<double> a0 = A00Inv*(b0 - A01*a1);

      ans.setBlock(0, a0);
      ans.setBlock(1, a1);

      cout << "x=" << x << std::endl;
      cout << "ans=" << ans << std::endl;

      VectorSpace<double> monoSpace 
        = type.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRows);
      PartitionedToMonolithicConverter converter(domain, isBCCol, monoSpace);

      Vector<double> monoVec = monoSpace.createMember();
      converter.convert(x, monoVec);
      
      cout << "monoVec=" << monoVec << std::endl;

      double err = (x-ans).normInf();
      double err0 = (x.getBlock(0)-ans.getBlock(0)).normInf();
      double err1 = (x.getBlock(1)-ans.getBlock(1)).normInf();

      cout << "error = " << err << std::endl;
      cout << "error0 = " << err0 << std::endl;
      cout << "error1 = " << err1 << std::endl;
      double tol = 1.0e-10;

      
      if (err > tol)
      {
        cout << "Poisson solve test FAILED" << std::endl;
      }
      else
      {
        cout << "Poisson solve test PASSED" << std::endl;
      }
#endif
    }

  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
  }
}

