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

#include "SundanceHermiteSpectralBasis.hpp"
#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMap.hpp"

using namespace std;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


HermiteSpectralBasis::HermiteSpectralBasis(int dim, int order)
  : SpectralBasisBase(),
    basis_(),
    dim_(dim),
    order_(order),
    maxterms_(-1),
    cijk_()
{
  Chaos Hermite(dim_, order_);
  maxterms_ = Hermite.tnterms();
  basis_.resize(maxterms_);
  for (int i=0; i<maxterms_; i++)
    {
      basis_[i] = i;
    }  
   cijk_ = rcp(new cijk(dim_, order_));
}

HermiteSpectralBasis::HermiteSpectralBasis(int dim, int order, int nterms)
  : SpectralBasisBase(),
    basis_(),
    dim_(dim),
    order_(order),
    maxterms_(nterms),
    cijk_()
{
  basis_.resize(maxterms_);
  for (int i=0; i<maxterms_; i++)
    {
      basis_[i] = i;
    }
  cijk_ = rcp(new cijk(dim_, order_));
}

HermiteSpectralBasis::HermiteSpectralBasis(int dim, int order, const Array<int>& basisarray)
  : SpectralBasisBase(),
    basis_(),
    dim_(dim),
    order_(order),
    maxterms_(-1),
    cijk_()
{
  maxterms_ = basisarray.size();
  basis_.resize(maxterms_);
  for (int i=0; i<maxterms_; i++)
    basis_[i] = basisarray[i];
  cijk_ = rcp(new cijk(dim_, order_));
}


int HermiteSpectralBasis::getDim() const 
{
  return dim_;
}

int HermiteSpectralBasis::getOrder() const
{
  return order_;
}

int HermiteSpectralBasis::nterms() const
{
  return maxterms_;
}

int HermiteSpectralBasis::getElement(int i) const
{
  return basis_[i];
}

double HermiteSpectralBasis::expectation(int i, int j, int k)
{
   return cijk_->expectation(i,j,k);
}


string HermiteSpectralBasis::toString() const
{
  return "HermiteSpectralBasis(" + Teuchos::toString(getDim())
    + ", " + Teuchos::toString(getOrder()) + ")";
}


bool HermiteSpectralBasis::lessThan(const SpectralBasisBase* other) const
{
  if (typeid(*this).before(typeid(*other))) return true;
  if (typeid(*other).before(typeid(*this))) return false;

  if (getDim() < other->getDim()) return true;
  if (other->getDim() < getDim()) return false;

  if (getOrder() < other->getOrder()) return true;
  if (other->getOrder() < getOrder()) return false;

  return false;
}
