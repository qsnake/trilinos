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
#include "TSFGlobalAnd.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFSerialVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFVectorTester.hpp"


#ifdef _MSC_VER
#include "winmath.h"
#endif

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


bool runTest(int nProc, int rank, const VectorType<double>& vecType)
{
  int n = 4;
  int seed = 12345;
  for (int i=0; i<rank; i++) seed = (seed * 371761) % 5476181;
  cout << "seed = " << seed << std::endl;
  srand48(seed);
   
  int dimension = nProc*n;
  int low = n*rank;
  std::vector<int> localRows(n);
  for (int i=0; i<n; i++)
  {
    localRows[i] = low + i;
  }
   
  VectorSpace<double> space = vecType.createSpace(dimension, n, 
    &(localRows[0]), MPIComm::world());
   
  VectorTester<double> tester(space, TestSpecifier<double>(true, 1.0e-13, 1.0e-10));
   
  bool allPass = tester.runAllTests();
  return allPass;
}

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
  {
    GlobalMPISession session(&argc, &argv);
    int nProc = session.getNProc();
    int rank = session.getRank();

    VectorType<double> type1 = new EpetraVectorType();
    VectorType<double> type2 = new SerialVectorType();
      
    bool allPass = true;

    allPass = runTest(nProc, rank, type1);

    if (rank==0)
    {
      allPass = runTest(1, rank, type2) && allPass;
    }

    allPass = globalAnd(allPass);
      
    if (!allPass) 
    {
	Out::root() << "detected a test that FAILED" << std::endl;
	stat = -1;
    }
    else
    {
	Out::root() << "all tests PASSED" << std::endl;
    }


  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << std::endl;
    stat = -1;
  }
  return stat;
}

