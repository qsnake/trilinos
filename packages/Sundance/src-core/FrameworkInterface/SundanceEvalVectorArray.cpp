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


#include "SundanceEvalVectorArray.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;

EvalVectorArray::EvalVectorArray()
  : Array<RCP<EvalVector> >()
{;}

void EvalVectorArray::copy(const RCP<EvalVectorArray>& other)
{
  resize(other->size());

  for (int i=0; i<other->size(); i++)
    {
      (*this)[i]->copy((*other)[i]);
    }
}

void EvalVectorArray::steal(const RCP<EvalVectorArray>& other)
{
  resize(other->size());

  for (int i=0; i<other->size(); i++)
    {
      (*this)[i] = (*other)[i];
    }
}

ostream& EvalVectorArray::print(std::ostream& os, 
                                const SparsitySuperset* derivs) const
{
  Tabs tab;
  TEST_FOR_EXCEPTION(derivs->numDerivs() != this->size(),
                     InternalError,
                     "mismatch between deriv set size=" << derivs->numDerivs()
                     << "and result vector size " << this->size()
                     << "in EvalVectorArray::print");

  int maxlen = 25;
  for (int i=0; i<derivs->numDerivs(); i++)
    {
      int s = derivs->deriv(i).toString().length();
      if (s > maxlen) maxlen = s;
    }
  
  for (int i=0; i<derivs->numDerivs(); i++)
    {
      os << tab;
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs->deriv(i).toString() << "\t\t";
      (*this)[i]->print(os);
      os << std::endl;
    }
  return os;
}
