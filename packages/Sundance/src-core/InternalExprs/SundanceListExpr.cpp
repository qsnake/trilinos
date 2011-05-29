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

#include "SundanceListExpr.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

static Time& appendToListTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("append to list"); 
  return *rtn;
}

ListExpr::ListExpr()
  : ExprBase(), elements_()
{;}

ListExpr::ListExpr(const Array<Expr>& elements)
  : ExprBase(), elements_(elements)
{;}

void ListExpr::append(const Expr& expr)
{
  TimeMonitor timer(appendToListTimer());
  elements_.append(expr);
}

Expr ListExpr::flatten() const 
{
  Expr rtn = new ListExpr();

  for (int i=0; i<this->size(); i++)
    {
      Expr e = element(i).flatten();
      for (int j=0; j<e.size(); j++)
        {
          rtn.append(e[j]);
        }
    }

  return rtn;
}

Expr ListExpr::join(const Expr& other) const 
{
  Expr rtn = new ListExpr(elements_);
  
  for (int i=0; i<other.size(); i++)
    {
      rtn.append(other[i]);
    }

  return rtn;
}

int ListExpr::size() const
{
  return elements_.size();
}

int ListExpr::totalSize() const 
{
  int rtn = 0;

  for (int i=0; i<this->size(); i++)
    {
      rtn += elements_[i].totalSize();
    }

  return rtn;
}

std::ostream& ListExpr::toText(std::ostream& os, bool paren) const
{
  os << "{";
  for (int i=0; i<elements_.size(); i++)
    {
      elements_[i].ptr()->toText(os, paren);
      if (i < elements_.size()-1) os << ", ";
    }
  os << "}";
  return os;
}


XMLObject ListExpr::toXML() const 
{
  XMLObject rtn("ListExpr");
  for (int i=0; i<elements_.length(); i++)
    {
      rtn.addChild(elements_[i].toXML());
    }
  return rtn;
}


