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

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosThyraAdapter.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "TSFPreconditioner.hpp"
#include "TSFPreconditionerFactory.hpp"
#include "TSFParameterListPreconditionerFactory.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "TSFAztecSolver.hpp"
#include "TSFMatrixLaplacian1D.hpp"
#include "TSFLinearSolverBuilder.hpp"
#include "TSFInverseOperatorDecl.hpp"
#include "SundancePathUtils.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "TSFLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#include "TSFInverseOperatorImpl.hpp"

#endif

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
  {
    GlobalMPISession session(&argc, &argv);

    MPIComm::world().synchronize();

    VectorType<double> type = new EpetraVectorType();


    ParameterXMLFileReader reader("belos-ml.xml");
    ParameterList solverParams = reader.getParameters();

    /* create the range space  */
    int nLocalRows = 20;

    MatrixLaplacian1D builder(nLocalRows, type);

    LinearOperator<double> A = builder.getOp();

    Vector<double> x = A.domain().createMember();
    int myRank = MPIComm::world().getRank();
    int nProcs = MPIComm::world().getNProc();
      

    Thyra::randomize(-ST::one(),+ST::one(),x.ptr().ptr());
#ifdef TRILINOS_6
    if (myRank==0) x[0] = 0.0;
    if (myRank==nProcs-1) x[nProcs * nLocalRows - 1] = 0.0;
    // need to fix operator[] routine
#else
    if (myRank==0) x.setElement(0, 0);
    if (myRank==nProcs-1) x.setElement(nProcs * nLocalRows - 1, 0.0);
#endif
    cout << "input is " << std::endl;
    x.print(cout);
    Vector<double> b = A*x;

    cout << "rhs is " << std::endl;
    b.print(cout);




    LinearSolver<double> solver = LinearSolverBuilder::createSolver(solverParams);
    LinearOperator<double> AInv = inverse(A, solver);
    Vector<double> ans = AInv * b;

    cout << "answer is " << std::endl;
    ans.print(cout);
      
    double err = (x-ans).norm2()/((double) nProcs * nLocalRows);
    cout << "error norm = " << err << std::endl;

    if (err <= 1.0e-8)
    {
      cout << "Belos poisson solve test PASSED" << std::endl;
      return 0;
    }
    else
    {
      cout << "Belos poisson solve test FAILED" << std::endl;
      return 1;
    }
  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

