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

#include "SundanceTestFunctionStub.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



TestFunctionStub::TestFunctionStub(const std::string& name, 
  int tensorOrder,
  int dim, 
  const RCP<const TestFuncDataStub>& data)
  : SymbolicFunc(makeFuncID(tensorOrder), 
    rcp_dynamic_cast<const CommonFuncDataStub>(data)), data_(data)
{
  FunctionIdentifier myFid = fid();
  if (tensorOrder==0)
  {
    append(new TestFuncElement(data, name, "", myFid));
  }
  else if (tensorOrder==1)
  {
    for (int d=0; d<dim; d++)
    {
      std::string suffix="[" + Teuchos::toString(d) + "]";
      FunctionIdentifier fid = myFid.createComponent(d);
      append(new TestFuncElement(data, name, suffix, fid));
    }
  }
  else 
  {
    TEST_FOR_EXCEPTION(true, RuntimeError, "tensor order = " << tensorOrder
      << " not supported");
  }
}

TestFunctionStub::TestFunctionStub(const std::string& name, 
  const SpectralBasis& sbasis, int tensorOrder, int dim,
  const RCP<const TestFuncDataStub>& data)
  :  SymbolicFunc(makeFuncID(tensorOrder), 
    rcp_dynamic_cast<const CommonFuncDataStub>(data)), data_(data)
{
  Array<FunctionIdentifier> cFid(sbasis.nterms());

  for (int n=0; n<sbasis.nterms(); n++)
  {
    cFid[n] = makeFuncID(tensorOrder);
  }
  
  if (tensorOrder==0 || dim==1)
  {
    Array<Expr> coeffs(sbasis.nterms());
    for (int n=0; n<sbasis.nterms(); n++)
    {
      std::string suffix="";
      if (sbasis.nterms()>1) suffix = "[" + Teuchos::toString(n) + "]";
      coeffs[n] = new TestFuncElement(data, name, suffix, cFid[n]);
    }
    append(new SpectralExpr(sbasis, coeffs));
  }
  else if (tensorOrder==1)
  {
    for (int d=0; d<dim; d++)
    {
      std::string suffix="[" + Teuchos::toString(d) + "]";
      Array<Expr> coeffs(sbasis.nterms());
      for (int n=0; n<sbasis.nterms(); n++)
      {
        FunctionIdentifier fid = cFid[n].createComponent(d);
        if (sbasis.nterms()>1) suffix += "[" + Teuchos::toString(n) + "]";
        coeffs[n]= new TestFuncElement(data, name, suffix, fid);
      }
      append(new SpectralExpr(sbasis, coeffs));
    }
  }
  else 
  {
    TEST_FOR_EXCEPTION(true, RuntimeError, "tensor order = " << tensorOrder
      << " not supported");
  }
}




