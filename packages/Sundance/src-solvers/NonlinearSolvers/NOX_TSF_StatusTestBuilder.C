// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_TSF_StatusTestBuilder.H"         
#include "NOX_StatusTest_NormF.H"         
#include "NOX_StatusTest_NormUpdate.H"         
#include "NOX_StatusTest_SafeCombo.H"         
#include "NOX_StatusTest_MaxIters.H"         
#include "Teuchos_TestForException.hpp"   

using namespace NOX;
using namespace NOX::TSF;
using namespace Teuchos;

RCP<StatusTest::Generic> 
StatusTestBuilder::makeStatusTest(const ParameterList& params)
{
  TEST_FOR_EXCEPTION(!params.isSublist("Status Test"), runtime_error,
                     "did not find Status Test sublist in " << params);

  ParameterList testSublist = params.sublist("Status Test");

  double fTol = 1.0e-15;
  double dxTol = 1.0e-15;
  int maxiters = 20;
  if (testSublist.isParameter("Tolerance"))
    {
      fTol = getParameter<double>(testSublist, "Tolerance");
    }
  if (testSublist.isParameter("Residual Tolerance"))
    {
      fTol = getParameter<double>(testSublist, "Residual Tolerance");
    }
  if (testSublist.isParameter("Step Tolerance"))
    {
      dxTol = getParameter<double>(testSublist, "Step Tolerance");
    }
  if (testSublist.isParameter("Max Iterations"))
    {
      maxiters = getParameter<int>(testSublist, "Max Iterations");
    }

  RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(fTol));
  RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(maxiters));
  RCP<StatusTest::Generic> C = rcp(new StatusTest::NormUpdate(dxTol));
  RCP<StatusTest::Generic> AB 
    = rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
  RCP<StatusTest::Generic> ABC 
    = rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, AB, C));
  
  return ABC;
}



