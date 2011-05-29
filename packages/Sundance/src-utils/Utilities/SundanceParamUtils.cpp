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

#include "SundanceParamUtils.hpp"

using Teuchos::Array;
using Teuchos::ParameterList;

using std::ifstream;

namespace Sundance
{
ParameterList mergeParamLists(const ParameterList& pDef, 
  const ParameterList& pIn)
{
  ParameterList rtn = pDef;
  using namespace Teuchos;
  
  /* replace any defaults with overriden values */
  ParameterList::ConstIterator i;

  for (i=pDef.begin(); i!=pDef.end(); i++)
  {
    const ParameterEntry& eDef = pDef.entry(i);

    const std::string& name = pDef.name(i);
    const ParameterEntry* eIn = pIn.getEntryPtr(name);
    if (eIn != NULL)
    {
      if (eIn->isList() && eDef.isList())
      {
        ParameterList sub = mergeParamLists(
          getValue<ParameterList>(eDef),
          getValue<ParameterList>(*eIn));
        rtn.set(name, sub);

      }
      else if (eIn->isType<int>() && eDef.isType<int>())
      {
        rtn.set(name, Teuchos::any_cast<int>(eIn->getAny()));
      }
      else
      {
        TEST_FOR_EXCEPTION(eIn->isList() && !eDef.isList(), 
          std::runtime_error, "mismatched parameters in mergeParams()");
        TEST_FOR_EXCEPTION(!eIn->isList() && eDef.isList(), 
          std::runtime_error, "mismatched parameters in mergeParams()");
        TEST_FOR_EXCEPT(1);
      }
    }
    else
    {
    }
  }
  return rtn;
}
}
