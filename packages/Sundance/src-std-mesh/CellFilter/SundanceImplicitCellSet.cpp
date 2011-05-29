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

#include "SundanceImplicitCellSet.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

ImplicitCellSet::ImplicitCellSet(const Mesh& mesh, int cellDim,
                                 const CellType& cellType)
  : CellSetBase(mesh, cellDim, cellType),
    maxLID_(mesh.numCells(cellDim))
{;}

CellIterator ImplicitCellSet::begin() const
{
  return CellIterator(mesh(), dimension(), CellIterator::Begin);
}

CellIterator ImplicitCellSet::end() const
{
  return CellIterator(mesh(), dimension(), CellIterator::End);
}

void ImplicitCellSet::print(std::ostream& os) const 
{
  os << "ImplicitCellSet[dim=" << dimension() << ", type=" << cellType() << "]";
}

bool ImplicitCellSet::internalLessThan(const CellSetBase* other) const
{
  const ImplicitCellSet* e = dynamic_cast<const ImplicitCellSet*>(other);

  if (e==0) return false;

  bool rtn = maxLID_ < e->maxLID_;
  return rtn;
}
