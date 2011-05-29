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

#include "SundanceStokhosBasisWrapper.hpp"

#ifdef HAVE_SUNDANCE_STOKHOS

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMap.hpp"
#include "SundanceExceptions.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"


namespace Sundance
{
using namespace Teuchos;


StokhosBasisWrapper::StokhosBasisWrapper(
  const RCP<const PolyBasis>& basis)
  : basis_(basis), cijk_()
{
  fillCijk();
}

StokhosBasisWrapper::StokhosBasisWrapper(
  const RCP<const PolyBasis1D>& basis)
  : basis_(), cijk_()
{
  Array<RCP<const PolyBasis1D> > bases(1); 
  bases[0] = basis;
  basis_ = rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  fillCijk();
}

void StokhosBasisWrapper::fillCijk()
{
  cijk_ = basis_->computeTripleProductTensor(basis_->order());
}


double StokhosBasisWrapper::expectation(int i, int j, int k)
{
  int I;
  int J;
  int N = cijk_->num_values(k);
  for (int L=0; L<N; L++)
  {
    double v;
    cijk_->value(k, L, I, J, v);
    if (i==I && j==J) return v;
  }
  return 0.0;
}




bool StokhosBasisWrapper::lessThan(const SpectralBasisBase* other) const
{
  if (typeid(*this).before(typeid(*other))) return true;
  if (typeid(*other).before(typeid(*this))) return false;

  if (getDim() < other->getDim()) return true;
  if (other->getDim() < getDim()) return false;

  if (getOrder() < other->getOrder()) return true;
  if (other->getOrder() < getOrder()) return false;

  const PolyBasis* me = basis_.ptr().get();
  /* because the type comparisons have been neutral, this cast should not
   * fail */
  const StokhosBasisWrapper* otherGuy 
    = dynamic_cast<const StokhosBasisWrapper*>(other);
  TEST_FOR_EXCEPTION(otherGuy == 0, RuntimeError,
    "unexpected cast failure");

  const PolyBasis* you = otherGuy->basis_.ptr().get();
  if (typeid(*me).before(typeid(*you))) return true;

  return false;
}


}

#endif
