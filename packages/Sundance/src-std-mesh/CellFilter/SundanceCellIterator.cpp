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

#include "SundanceCellIterator.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

CellIterator::CellIterator()
  :  isImplicit_(true),
    currentLID_(-1),
    reorderer_(0),
     iter_(dummy().begin())
{;}

CellIterator::CellIterator(const CellIterator& other)
  :  isImplicit_(other.isImplicit_),
     currentLID_(other.currentLID_),
     reorderer_(other.reorderer_),
     iter_()
{
  if (!isImplicit_) iter_ = other.iter_;
}

CellIterator::CellIterator(const Mesh& mesh, 
                           int cellDim, 
                           CellIteratorPos pos)
  : isImplicit_(true),
    currentLID_(-1),
    reorderer_(mesh.reorderer()),
    iter_(dummy().begin())
{
  if (cellDim == mesh.spatialDim() && reorderer_ != 0)
    {
      switch(pos)
        {
        case Begin:
          currentLID_ = reorderer_->begin(); 
          break;
        case End:
          currentLID_ = mesh.numCells(cellDim);
        }
      SUNDANCE_OUT(mesh.verb() > 2, 
                   "created implicit cell iterator with LID=" << currentLID_);
    }
  else 
    {
      switch(pos)
        {
        case Begin:
          currentLID_ = 0;
          break;
        case End:
          currentLID_ = mesh.numCells(cellDim);
        }
      SUNDANCE_OUT(mesh.verb() > 2, 
                   "created implicit cell iterator with LID=" << currentLID_);
    }


}



CellIterator::CellIterator(const Set<int>* cells, CellIteratorPos pos)
  : isImplicit_(false),
    currentLID_(-1),
    reorderer_(0),
    iter_(dummy().begin())
{
  switch(pos)
  {
    case Begin:
      iter_ = cells->begin();
      break;
    case End:
      iter_ = cells->end();
      break;
    default:
      TEST_FOR_EXCEPT(1);
  }
}



CellIterator& CellIterator::operator=(const CellIterator& other)
{
  if (*this!=other) 
  {
    isImplicit_ = other.isImplicit_;
    currentLID_ = other.currentLID_;
    reorderer_=other.reorderer_;
    if (!isImplicit_) iter_ = other.iter_;
  }
  return *this;
}

    
