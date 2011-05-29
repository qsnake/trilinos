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
#include "TSFGlobalAnd.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorOpsDecl.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFSerialVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFMultiVectorOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFVectorImpl.hpp"
#include "TSFVectorOpsImpl.hpp"
#include "TSFMultiVectorOperatorImpl.hpp"
#endif


using namespace TSFExtended;
using namespace TSFExtendedOps;
using Thyra::TestSpecifier;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
    {
      GlobalMPISession session(&argc, &argv);
 
      /* create a distributed vector space for the multivector's vectors */
      VectorType<double> rowType = new EpetraVectorType();
      int nLocalRows = 10;
      VectorSpace<double> space 
        = rowType.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRows);
      
      /* create a replicated vector space for the small space of columns */
      int nVecs = 3;
      VectorType<double> colType = new SerialVectorType();
      VectorSpace<double> replSpace 
        = colType.createEvenlyPartitionedSpace(MPIComm::world(), nVecs);

      /* create some random vectors */
      Teuchos::Array<Vector<double> > vecs(nVecs);
      for (int i=0; i<nVecs; i++)
      {
        vecs[i] = space.createMember();
        /* do the operation elementwise */
        randomize(vecs[i]);
      }

      /* Test multiplication by a multivector operator. We will compute
       * y1 by directly summing columns, and y2 by applying the operator */
      LinearOperator<double> A = multiVectorOperator<double>(vecs, replSpace);

      
      Vector<double> y1 = space.createMember();
      Vector<double> y2 = space.createMember();
      y1.zero(); 
      y2.zero(); 
      

      /* Sum columns, putting the weights into x */
      Vector<double> x = replSpace.createMember();
      SequentialIterator<double> j;
      for (j=replSpace.begin(); j != replSpace.end(); j++)
        {
          x[j] = 2.0*(drand48()-0.5);
          y1 = y1 + x[j] * vecs[j.globalIndex()];
        }  
      
      /* Apply the operator to the vector of weights */
      y2 = A * x;
      
      double errA = (y1-y2).normInf();

      Out::root() << "error in A*x = " << errA << std::endl;


      /* Now test z = A^T * y */
      LinearOperator<double> At = A.transpose();
      
      Vector<double> z1 = replSpace.createMember();
      z1.zero();
      Vector<double> z2 = replSpace.createMember();
      z2.zero();

      Vector<double> y = y1.copy();

      /* compute by vectorwise multiplication */
      for (j=replSpace.begin(); j != replSpace.end(); j++)
      {
        z1[j] = vecs[j.globalIndex()].dot(y);
      }
      /* compute with operator */
      z2 = At * y;
      

      double errAt = (z1-z2).normInf();
      Out::root() << "error in At*y = " << errA << std::endl;

      double tol = 1.0e-13;
      bool pass = errA + errAt < tol;
      pass = globalAnd(pass);
      if (pass)
        {
          Out::root() << "multivector op test PASSED" << std::endl;
        }
      else
        {
          stat = -1;
          Out::root() << "multivector op test FAILED" << std::endl;
        }
    }
  catch(std::exception& e)
    {
      stat = -1;
      std::cerr << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



