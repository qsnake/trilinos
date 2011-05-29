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


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(Sundance::searchForFile("SolverParameters/userPrecParams.xml"));
#else
      ParameterXMLFileReader reader("userPrecParams.xml");
#endif

    ParameterList solverParams = reader.getParameters();
    ParameterList innerSolverParams = solverParams.sublist("Inner Solve");
    ParameterList outerSolverParams = solverParams.sublist("Outer Solve");

    /* create the range space  */
    int nLocalRows = solverParams.get<int>("nLocal");

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

    cout << "x=" << std::endl;
    x.print(cout);
      
    Vector<double> y = A*x;
    cout << "y=" << std::endl;
    y.print(cout);

    Vector<double> ans = A.range().createMember();

    LinearSolver<double> innerSolver 
      = LinearSolverBuilder::createSolver(innerSolverParams);

    LinearSolver<double> outerSolver 
      = LinearSolverBuilder::createSolver(outerSolverParams);

    /* call the setUserPrec() function to set the operator and solver 
     * to be used for preconditioning */
    outerSolver.setUserPrec(A, innerSolver);

    LinearOperator<double> AInv = inverse(A, outerSolver);
      

    ans = AInv * y;

    //      SolverState<double> state = solver.solve(A, y, ans);
      

      
    //      cout << state << std::endl;

    cout << "answer is " << std::endl;
    ans.print(cout);
      
    double err = (x-ans).norm2();
    cout << "error norm = " << err << std::endl;

    double tol = 1.0e-8;
    if (err > tol)
    {
      cout << "User-defined preconditioner test FAILED" << std::endl;
      return 1;
    }
    else
    {
      cout << "User-defined preconditioner test PASSED" << std::endl;
      return 0;
    }
  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

