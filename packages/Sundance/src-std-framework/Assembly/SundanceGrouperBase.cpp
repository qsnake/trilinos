/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceGrouperBase.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceUnknownParameterElement.hpp"
#include "SundanceTestFunction.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;



void GrouperBase::setVerbosity(
  int setupVerb,
  int integrationVerb,
  int transformVerb)
{
  setupVerb_ = setupVerb;
  integrationVerb_ = integrationVerb;
  transformVerb_ = transformVerb;
}



void GrouperBase::extractWeakForm(const EquationSet& eqn,
  const MultipleDeriv& functionalDeriv,
  BasisFamily& varBasis, 
  BasisFamily& unkBasis,
  MultiIndex& miVar, MultiIndex& miUnk,
  int& rawVarID, int& rawUnkID,  
  int& reducedVarID, int& reducedUnkID,  
  int& testBlock, int& unkBlock, 
  int& rawParamID, int& reducedParamID, 
  bool& isOneForm, bool& hasParam) const
{
  Tabs tab0(0);

  MultipleDeriv::const_iterator iter;

  isOneForm = false;  
  hasParam = false;

  if (functionalDeriv.size()==0) return;

  TEST_FOR_EXCEPTION(functionalDeriv.size() > 2, InternalError,
    "GrouperBase::extractWeakForm detected a functional "
    "derivative of order > 2: " 
    << functionalDeriv.toString());

  bool foundUnk = false;
  bool foundVar = false;

  SUNDANCE_MSG2(setupVerb(), 
    tab0 << "extracting weak form for functional derivative " 
    << functionalDeriv);


  for (iter = functionalDeriv.begin(); iter != functionalDeriv.end(); iter++)
  {
    Tabs tab;
    const Deriv& d = *iter;
      
    TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), InternalError,
      "GrouperBase::extractWeakForm "
      "detected a non-functional derivative: "
      << functionalDeriv.toString());
      
    const FunctionIdentifier& fid = d.fid();
      
    const SymbolicFuncElement* s = d.symbFuncElem();

    TEST_FOR_EXCEPTION(s==0, InternalError, 
      "GrouperBase::extractWeakForm failed to cast "
      "function to SymbolicFuncElement");
      

    int dofID = fid.dofID();
    int myIndex = fid.componentIndex();

    if (!foundVar && eqn.hasVarID(dofID))
    {
      TEST_FOR_EXCEPTION(d.isParameter(), InternalError,
        "Parameter not expected here");
      foundVar = true;
      reducedVarID = eqn.reducedVarID(dofID);
      rawVarID = dofID;
      testBlock = eqn.blockForVarID(dofID);

      SUNDANCE_MSG2(setupVerb(), 
        tab << "found varID=" << reducedVarID);

      const UnknownFuncElement* u
        = dynamic_cast<const UnknownFuncElement*>(s);

      const TestFuncElement* t
        = dynamic_cast<const TestFuncElement*>(s);

      TEST_FOR_EXCEPTION(u==0 && t==0, InternalError, 
        "GrouperBase::extractWeakForm could not cast "
        "variational function to either an "
        "UnknownFuncElement or a TestFuncElement");

      if (t != 0) 
      {
        varBasis = TestFunctionData::getData(t)->basis()[myIndex];
      }
      else 
      {
        varBasis = UnknownFunctionData::getData(u)->basis()[myIndex];
      }
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found varBasis=" << varBasis);

      miVar = d.opOnFunc().mi();
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found var multi index=" << miVar.toString());
    }
    else if (eqn.hasFixedParamID(dofID))
    {
      const UnknownParameterElement* upe
        = dynamic_cast<const UnknownParameterElement*>(s);
      TEST_FOR_EXCEPTION(upe==0, InternalError, 
        "GrouperBase::extractWeakForm could not cast "
        "unknown parameter to UnknownParameterElement");
      hasParam = true;
      rawParamID = dofID;
      reducedParamID = eqn.reducedFixedParamID(dofID);
    }
    else
    {
      TEST_FOR_EXCEPTION(d.isParameter(), InternalError,
        "Parameter not expected here");
      const UnknownFuncElement* u
        = dynamic_cast<const UnknownFuncElement*>(s);
      TEST_FOR_EXCEPTION(u==0, InternalError, 
        "GrouperBase::extractWeakForm could not cast "
        "unknown function to UnknownFuncElement");
      foundUnk = true;
      reducedUnkID = eqn.reducedUnkID(dofID);
      rawUnkID = dofID;
      unkBlock = eqn.blockForUnkID(dofID);

      SUNDANCE_MSG2(setupVerb(), 
        tab << "found reducedUnkID=" << reducedUnkID);

      unkBasis = UnknownFunctionData::getData(u)->basis()[myIndex];
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found unkBasis=" << unkBasis);

      miUnk = d.opOnFunc().mi();
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found unk multi index=" << miUnk.toString());
    }
  }

  if (!foundUnk) isOneForm = true;
}
