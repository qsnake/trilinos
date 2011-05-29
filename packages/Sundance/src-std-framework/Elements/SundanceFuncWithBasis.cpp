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

#include "SundanceFuncWithBasis.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"


namespace Sundance
{
using namespace Teuchos;

std::string describeFunction(const Expr& f)
{
  TEST_FOR_EXCEPT(f.ptr().get()==0);

  if (f.size() == 1)
  {
    const FuncElementBase* fe = dynamic_cast<const FuncElementBase*>(f[0].ptr().get());
    TEST_FOR_EXCEPTION(fe==0, RuntimeError, "expected a FuncElementBase, "
      "found " << typeid(*fe).name());
    
    const UnknownFuncElement* u = dynamic_cast<const UnknownFuncElement*>(f[0].ptr().get());

    const TestFuncElement* t = dynamic_cast<const TestFuncElement*>(f[0].ptr().get());

    const DiscreteFuncElement* d = dynamic_cast<const DiscreteFuncElement*>(f[0].ptr().get());

    std::string type;
    if (t != 0) 
    {
      type = "TFElem";
    }
    else if (u != 0) 
    {
      type = "UFElem";
    }
    else if (d != 0)
    {
      type = "DFElem";
    }
    else
    {
      TEST_FOR_EXCEPTION(true, RuntimeError, "unrecognized function " 
        << f[0]);
    }

    std::string rtn = type + "[name=" + fe->name() + ", fid=" + fe->fid().toString() + "]";
    return rtn;
      
  }
  else
  {
    std::string rtn = "{";
    for (int i=0; i<f.size(); i++)
    {
      if (i != 0) rtn += ", ";
      rtn += describeFunction(f[i]);
    }
    rtn += "}";
    return rtn;
  }
}
}
