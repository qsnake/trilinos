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
#include "TSFVectorDecl.hpp"
#include "TSFDefaultLTIProblemFactory.hpp"
#include "TSFCommonOperatorsDecl.hpp"
#include "TSFLinearCombinationDecl.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFRandomSparseMatrix.hpp"
#include "TSFProductVectorSpaceDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFCommonOperatorsDecl.hpp"
#endif


STREAM_OUT(Vector<double>)
//using namespace Teuchos;
  using namespace TSFExtended;
using namespace TSFExtendedOps;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
  {
    GlobalMPISession session(&argc, &argv);
 

    VectorType<double> type = new EpetraVectorType();

    int nLocalRows = 4;
      
    double onProcDensity = 0.5;
    double offProcDensity = 0.1;

    RandomSparseMatrix<double> ABuilder(nLocalRows, nLocalRows, 
      onProcDensity, offProcDensity, type);

    int nSteps = 3;

    /* Create the operator that advances a single timestep */
    LinearOperator<double> A = ABuilder.getOp();
    /* Create the operator that produces observables from a state.
     * For this test, we observe all states so C=I */
    LinearOperator<double> C = identityOperator<double>(A.domain());

    /* Create the factory for building the multi-timestep operators */
    RCP<DefaultLTIProblemFactory<double> > fact
      = rcp(new DefaultLTIProblemFactory<double>(nSteps));
    fact->init(A, C);


    /* --------------------------------------------------------------- 
     *                      Compute forward solve        
     *
     * We'll do this in two ways: first, by applying the bigAInv 
     * operator; second, by writing out the timestepping loop. 
     * These are completely equivalent, so the computations had
     * better give the same results!
     * --------------------------------------------------------------- */
      
    /* bigAInv is the n-step integration operator */
    LinearOperator<double> bigAInv = fact->getBigAInv();
    /* bigF pads a vector with zeros everywhere but the first block */
    LinearOperator<double> bigF = fact->getBigF();

    /* Create a random initial state */
    Vector<double> x0 = bigF.domain().createMember();
    Thyra::randomize(-1.0, 1.0, x0.ptr().get());
      
    /* Do the timestepping with a single operator application */
    Vector<double> bigX = bigF * x0;
    Vector<double> bigY = bigAInv * bigX;

    /* Extract the final timestep */
    Vector<double> xf = bigY.getBlock(nSteps-1);

    /* Now, repeat the timestepping with hand-coded loop */
    Vector<double> y;
    for (int i=0; i<nSteps; i++)
    {
      if (i==0)
      {
        y = x0.getBlock(0).copy();
      }
      else
      {
        y = A*y;
      }
      std::cout << "i=" << i << std::endl 
                << " y=" << y << std::endl
                << " bigY=" << bigY.getBlock(i) << std::endl;
    }

    /* Compare results of operator-notation and hand-coded procedures.
     * These should be the same.*/
    double errFwd = (y-xf).norm2();
      
    std::cout << "forward solve err = " << errFwd << std::endl;

    /* --------------------------------------------------------------- 
     *                      Compute adjoint solve        
     *
     * We'll do this in two ways: first, by applying the bigAInvT 
     * operator; second, by writing out the backwards 
     * timestepping loop. These are completely equivalent, so the 
     * computations had better give the same results!
     * --------------------------------------------------------------- */

    bigX.zero();
    xf = bigX.getBlock(nSteps-1);
    Thyra::randomize(-1.0, 1.0, xf.ptr().get());

    bigY = bigAInv.transpose() * bigX;

    for (int i=nSteps-1; i>=0; i--)
    {
      if (i==nSteps-1)
      {
        y = xf.copy();
      }
      else
      {
        y = A.transpose()*y;
      }
      std::cout << "i=" << i << std::endl 
                << " y=" << y << std::endl
                << " bigY=" << bigY.getBlock(i) << std::endl;
    }

    x0 = bigY.getBlock(0);
    double errAdj = (y-x0).norm2();

    std::cout << "adjoint solve err = " << errAdj << std::endl;


    /* ----------------------------------------------------------- 
     *
     * Apply the Hessian 
     *
     * ----------------------------------------------------------- */
    LinearOperator<double> H = fact->getH();
    Vector<double> u0 = bigF.domain().createMember();
    std::cerr << "bigF.domain()=" << bigF.domain().description() << std::endl;
    std::cerr << "bigF.range()=" << bigF.range().description() << std::endl;


    Thyra::randomize(-1.0, 1.0, x0.ptr().get());

    Vector<double> z = H*u0;
      
    double tol = 1.0e-13;
    if (std::max(errAdj, errFwd) < tol)
    {
      std::cerr << "LTI test PASSED" << std::endl;
    }
    else
    {
      stat = -1;
      std::cerr << "LTI test FAILED" << std::endl;
    }
  }
  catch(std::exception& e)
  {
    stat = -1;
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
  return stat;
}



