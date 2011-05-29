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

#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



DiscreteFunctionStub::DiscreteFunctionStub(const std::string& name, 
  int tensorOrder,
  int dim, 
  const RCP<DiscreteFuncDataStub>& data,
  int listIndex)
  : ListExpr(), data_(data)
{
  initTensor(name, tensorOrder, dim, data,  listIndex);
}

void DiscreteFunctionStub::initTensor(const std::string& name, 
  int tensorOrder,
  int dim, 
  const RCP<DiscreteFuncDataStub>& data,
  int listIndex)
{
  FunctionIdentifier myFid = makeFuncID(tensorOrder);
  if (tensorOrder==0)
  {
    append(new DiscreteFuncElement(data, name, "", myFid, listIndex));
  }
  else if (tensorOrder==1)
  {
    for (int d=0; d<dim; d++)
    {
      std::string suffix="[" + Teuchos::toString(d) + "]";
      FunctionIdentifier fid = myFid.createComponent(d);
      append(new DiscreteFuncElement(data, name, suffix, fid, listIndex));
    }
  }
  else 
  {
    TEST_FOR_EXCEPTION(true, RuntimeError, "tensor order = " << tensorOrder
      << " not supported");
  }
}


void DiscreteFunctionStub::initTensorSpectral(const std::string& name, 
  const SpectralBasis& sbasis, 
  int tensorOrder,
  int dim, 
  const RCP<DiscreteFuncDataStub>& data,
  int listIndexOffset)
{
  Array<FunctionIdentifier> cFid(sbasis.nterms());
  Array<int> listIndex(sbasis.nterms());

  for (int n=0; n<sbasis.nterms(); n++)
  {
    cFid[n] = makeFuncID(tensorOrder);
    listIndex[n] = listIndexOffset + n;
  }
  
  if (tensorOrder==0 || dim==1)
  {
    Array<Expr> coeffs(sbasis.nterms());
    for (int n=0; n<sbasis.nterms(); n++)
    {
      std::string suffix="";
      if (sbasis.nterms()>1) suffix = "[" + Teuchos::toString(n) + "]";
      coeffs[n] = new DiscreteFuncElement(data, name, suffix, cFid[n], listIndex[n]);
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
        coeffs[n]= new DiscreteFuncElement(data, name, suffix, fid, listIndex[n]);
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




DiscreteFunctionStub::DiscreteFunctionStub(const Array<string>& name, 
  const Array<std::pair<int,int> >& tensorStructure,
  const RCP<DiscreteFuncDataStub>& data)
  : ListExpr(), data_(data)
{
  TEST_FOR_EXCEPT(name.size() != tensorStructure.size() && name.size()!=1);

  if (tensorStructure.size()==1)
  {
    int tensorOrder = tensorStructure[0].first;
    int dim = tensorStructure[0].second;
    initTensor(name[0], tensorOrder, dim, data, 0);
  }
  else
  {
    for (int i=0; i<tensorStructure.size(); i++)
    {
      std::string nm;
      if (name.size()==1) nm = name[0] + "[" + Teuchos::toString(i) + "]";
      else nm = name[i];
      append(new DiscreteFunctionStub(
               nm,
               tensorStructure[i].first,
               tensorStructure[i].second,
               data, i));
    }
  }
}



DiscreteFunctionStub::DiscreteFunctionStub(const std::string& name, 
  const SpectralBasis& sbasis, int tensorOrder, int dim,
  const RCP<DiscreteFuncDataStub>& data,
  int listIndex)
  : ListExpr(), data_(data)
{
  initTensorSpectral(name, sbasis, tensorOrder, dim, data, listIndex);
}



DiscreteFunctionStub::DiscreteFunctionStub(const Array<string>& name, 
  const SpectralBasis& sbasis,  
  const Array<std::pair<int,int> >& tensorStructure,
  const RCP<DiscreteFuncDataStub>& data)
  : ListExpr(), data_(data)
{
  TEST_FOR_EXCEPT(name.size() != tensorStructure.size());
   if (tensorStructure.size()==1)
  {
    int tensorOrder = tensorStructure[0].first;
    int dim = tensorStructure[0].second;
    initTensorSpectral(name[0], sbasis, tensorOrder, dim, data, 0);
  }
  else
  {
    for (int i=0; i<tensorStructure.size(); i++)
    {
      std::string nm;
      if (name.size()==1) nm = name[0] + "[" + Teuchos::toString(i) + "]";
      else nm = name[i];
      append(new DiscreteFunctionStub(
               nm,
               sbasis,
               tensorStructure[i].first,
               tensorStructure[i].second,
               data, i*sbasis.nterms()));
    }
  }
}


