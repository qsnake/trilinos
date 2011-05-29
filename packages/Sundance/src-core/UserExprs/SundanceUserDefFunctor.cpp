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

#include "SundanceUserDefFunctor.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


UserDefFunctor::UserDefFunctor(const std::string& name, 
                               int domainDim, 
                               int rangeDim)
  :  name_(name), 
     elemNames_(), 
     domainDim_(domainDim), 
     rangeDim_(rangeDim)
{
  TEST_FOR_EXCEPT(domainDim_ <= 0);
  TEST_FOR_EXCEPT(rangeDim_ <= 0);

  elemNames_.resize(rangeDim_);

  if (rangeDim_==1) 
    {
      elemNames_[0] = name_;
    }
  else
    {
      for (int i=0; i<rangeDim_; i++)
        {
          elemNames_[i] = name_ + "[" + Teuchos::toString(i) + "]";
        }
    }
}





void UserDefFunctor::eval0(const Array<double>& in, double* out) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "eval0() function not supported for functor " << name_);
}

void UserDefFunctor::eval1(const Array<double>& in, double* out, double* outDerivs) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "eval1() function not supported for functor " << name_);
}

