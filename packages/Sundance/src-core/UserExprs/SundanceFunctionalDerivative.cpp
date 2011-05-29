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

#include "SundanceFunctionalDerivative.hpp"

#include "SundanceSymbolicTransformation.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceExplicitFunctionalDerivativeElement.hpp"

using namespace Sundance;
using namespace Teuchos;

using namespace Sundance;

Expr Sundance::FunctionalDerivative(const Expr& F, const Expr& u)
{
  if (F.size() == 1 && u.size()==1)
  {
    RCP<ScalarExpr> arg = SymbolicTransformation::getScalar(F[0]);

    const UnknownFuncElement* uPtr
      = dynamic_cast<const UnknownFuncElement*>(u[0].ptr().get());
    
    Deriv fd = funcDeriv(uPtr);
    
    return new ExplicitFunctionalDerivativeElement(arg, fd);
  }
  else if (F.size() > 1)
  {
    Array<Expr> rtnList(F.size());
    for (int i=0; i<rtnList.size(); i++)
    {
      rtnList[i] = FunctionalDerivative(F[i], u);
    }
    return new ListExpr(rtnList);
  }
  else 
  {
    Array<Expr> rtnList(u.size());
    for (int i=0; i<rtnList.size(); i++)
    {
      rtnList[i] = FunctionalDerivative(F, u[i]);
    }
    return new ListExpr(rtnList);
  }
  
}

