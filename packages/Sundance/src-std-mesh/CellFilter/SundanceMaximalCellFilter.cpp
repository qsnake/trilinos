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

#include "SundanceMaximalCellFilter.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

MaximalCellFilter::MaximalCellFilter()
  : CellFilterBase()
{
  setName("MaximalCells");
}

CellSet MaximalCellFilter::internalGetCells(const Mesh& mesh) const
{
  return new ImplicitCellSet(mesh, mesh.spatialDim(), 
                             mesh.cellType(mesh.spatialDim()));
}

int MaximalCellFilter::dimension(const Mesh& mesh) const 
{
  return mesh.spatialDim();
}


bool MaximalCellFilter::lessThan(const CellFilterStub* /* other */) const
{
  /* Maximal cell sets always win a comparison */
  return false;
}



