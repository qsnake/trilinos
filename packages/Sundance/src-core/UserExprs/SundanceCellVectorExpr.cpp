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

#include "SundanceCellVectorExpr.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



namespace Sundance
{

Expr CellNormalExpr(int dimension, const std::string& name)
{
  Array<Expr> comps(dimension);
  for (int i=0; i<dimension; i++)
  {
    comps[i] = new CellVectorExpr(i, dimension, name + "["
      + Teuchos::toString(i) + "]");
  }
  return new ListExpr(comps);
}


Expr CellTangentExpr(int dimension, const std::string& name)
{
  Array<Expr> comp(dimension);
  for (int i=0; i<dimension; i++)
  {
    comp[i] = new CellVectorExpr(0, i, dimension, name + "("
        + Teuchos::toString(i) + ")");
  }
  return new ListExpr(comp);
}

}

CellVectorExpr::CellVectorExpr(int tangentBasisIndex, 
			       int tangentComponentIndex,
			       int dim,
			       const std::string& name)
  : EvaluatableExpr(), name_(name), dim_(dim), type_(CellTangentSpace),
    basisMemberIndex_(tangentBasisIndex), 
    componentIndex_(tangentComponentIndex)
{}


CellVectorExpr::CellVectorExpr(int normalComponentIndex, int dim,
  const std::string& name)
  : EvaluatableExpr(), name_(name), dim_(dim), type_(CellNormalVector),
    basisMemberIndex_(-1), componentIndex_(normalComponentIndex)
{}

bool CellVectorExpr::lessThan(const ScalarExpr* other) const
{
  const CellVectorExpr* f = dynamic_cast<const CellVectorExpr*>(other);
  TEST_FOR_EXCEPTION(f==0, InternalError, "cast should never fail at this point");
  if (type_ < f->type_) return true;
  if (type_ > f->type_) return false;
  if (dim_ < f->dim_) return true;
  if (dim_ > f->dim_) return false;
  if (basisMemberIndex_ < f->basisMemberIndex_) return true;
  if (basisMemberIndex_ > f->basisMemberIndex_) return false;
  return componentIndex_ < f->componentIndex_;
}


XMLObject CellVectorExpr::toXML() const 
{
  XMLObject rtn("CellVectorExpr");
  rtn.addAttribute("name", name_);
  return rtn;
}



Set<MultipleDeriv> 
CellVectorExpr::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  
  if (order==0) rtn.put(MultipleDeriv());
  
  return rtn;
}




std::ostream& CellVectorExpr::toText(std::ostream& os, bool paren) const
{
  os << name();
  return os;
}




