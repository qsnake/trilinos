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

#include "SundanceSymbolicFunc.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceZeroExpr.hpp"

#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

SymbolicFunc::SymbolicFunc(const FunctionWithID& fid, 
    const RCP<const CommonFuncDataStub>& data)
  : ListExpr(), FunctionWithID(fid), commonData_(data)
{}


void SymbolicFunc::substituteZero() const 
{
  for (int i=0; i<this->size(); i++)
    {
      const SymbolicFuncElement* u 
        = dynamic_cast<const SymbolicFuncElement*>(element(i).ptr().get());
      TEST_FOR_EXCEPTION(u==0, InternalError, 
                         "Non-symbolic function "
                         << element(i).toString() 
                         << " detected in SymbolicFunc::substituteZero()");
      u->substituteZero();
    }
}

void SymbolicFunc
::substituteFunction(const RCP<DiscreteFunctionStub>& u0) const
{
  TEST_FOR_EXCEPTION(this->size() != u0->size(), InternalError,
                     "Mismatch between sizes of symbolic " << toString()
                     << " and discrete func " << u0->toString()
                     << " in substituteFunction()");

  for (int i=0; i<this->size(); i++)
    {
      const SymbolicFuncElement* u 
        = dynamic_cast<const SymbolicFuncElement*>(element(i).ptr().get());
      TEST_FOR_EXCEPTION(u==0, InternalError, 
                         "Non-symbolic function "
                         << element(i).toString() 
                         << " detected in SymbolicFunc::substituteFunction()");

      RCP<DiscreteFuncElement> df 
        = rcp_dynamic_cast<DiscreteFuncElement>(u0->element(i).ptr());
      u->substituteFunction(df);
    }
}

