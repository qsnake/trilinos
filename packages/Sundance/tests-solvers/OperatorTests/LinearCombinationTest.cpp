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
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearCombinationTester.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#endif



using namespace Teuchos;
using namespace std;
using namespace TSFExtended;
using namespace TSFExtendedOps;
using Thyra::TestSpecifier;

int main(int argc, char *argv[]) 
{
  int stat = 0;

  try
    {
      GlobalMPISession session(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      int nLocalRows = 10;
      
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;

      LinearCombinationTester<double> tester(nLocalRows,
                                             onProcDensity,
                                             offProcDensity,
                                             type,
                                             TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      
      bool allPass = tester.runAllTests();
      if (!allPass) stat = -1;
    }
  catch(std::exception& e)
    {
      stat = -1;
      std::cerr << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



