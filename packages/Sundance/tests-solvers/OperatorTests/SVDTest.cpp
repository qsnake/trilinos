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
#include "TSFSerialVectorType.hpp"
#include "TSFDenseSerialMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorOpsDecl.hpp"
#include "TSFSimpleDiagonalOpDecl.hpp"
#include "TSFSimpleComposedOpDecl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "SundanceOut.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFVectorOpsImpl.hpp"
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimpleDiagonalOpImpl.hpp"
#include "TSFSimpleComposedOpImpl.hpp"
#endif



using namespace Teuchos;
using namespace Sundance;
using namespace TSFExtended;
using namespace TSFExtendedOps;
using std::endl;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
    {
      GlobalMPISession session(&argc, &argv);

      VectorType<double> vecType = new SerialVectorType();

      VectorSpace<double> domain 
        = vecType.createEvenlyPartitionedSpace(MPIComm::world(), 3);
      VectorSpace<double> range
        = vecType.createEvenlyPartitionedSpace(MPIComm::world(), 5);

      RCP<MatrixFactory<double> > mf 
        = vecType.createMatrixFactory(domain, range);

      LinearOperator<double> A = mf->createMatrix();
      RCP<DenseSerialMatrix> APtr 
        = rcp_dynamic_cast<DenseSerialMatrix>(A.ptr());

      APtr->setRow(0, tuple(1.0,   2.0,  3.0));
      APtr->setRow(1, tuple(4.0,   5.0,  6.0));
      APtr->setRow(2, tuple(7.0,   8.0,  9.0));
      APtr->setRow(3, tuple(10.0, 11.0,  12.0));
      APtr->setRow(4, tuple(13.0, 14.0,  15.0));

      Out::os() << "A = " << std::endl;
      A.setVerbosity(10);
      Out::os() << A << std::endl;

      LinearOperator<double> U;
      LinearOperator<double> Vt;
      Vector<double> sigma;

      denseSVD(A, U, sigma, Vt);

      Out::os() << "U = " << std::endl;
      U.setVerbosity(10);
      Out::os() << U << std::endl;

      Out::os() << "sigma = " << std::endl;
      Out::os() << sigma << std::endl; 

      Out::os() << "Vt = " << std::endl;
      Vt.setVerbosity(10);
      Out::os() << Vt << std::endl;

      int nSamples = 10;
      bool allOK = true;
      double tol = 1.0e-13;
      for (int i=0; i<nSamples; i++)
      {
        Out::os() << "Sample #" << i << " of " << nSamples << std::endl;
        Vector<double> x = domain.createMember();
        randomize(x);
        
        U.setVerbosity(0);
        Vt.setVerbosity(0);
        A.setVerbosity(0);
        
        LinearOperator<double> Sigma = diagonalOperator(sigma);
        
        Vector<double> z = (U * Sigma * Vt)*x - A*x;
        double ez = z.norm2();
        Out::os() << "|| (U Sigma Vt - A)*x || = " << ez << std::endl;
        
        Vector<double> y = (U.transpose() * U)*x - x;
        double ey = y.norm2();
        Out::os() << "|| (U^T U - I)*x || = " << ey << std::endl;
        
        Vector<double> w = (Vt * Vt.transpose())*x - x;
        double ew = w.norm2();
        Out::os() << "|| (V^T*V - I)*x || = " << ew << std::endl;
        if (ew > tol || ez > tol || ey > tol) allOK = false;
      }

      if (allOK)
      {
        Out::os() << "SVD test PASSED" << std::endl;
      }
      else
      {
        Out::os() << "SVD test FAILED" << std::endl;
        stat = -1;
      }

      
    }
  catch(std::exception& e)
    {
      stat = -1;
      Out::os() << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



