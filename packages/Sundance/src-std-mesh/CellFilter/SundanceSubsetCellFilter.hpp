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

#ifndef SUNDANCE_SUBSETCELLFILTER_H
#define SUNDANCE_SUBSETCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellFilterBase.hpp"
#include "SundanceCellPredicate.hpp"

namespace Sundance
{
using namespace Teuchos;

/** */
class SubsetCellFilter : public CellFilterBase 
{
public:
  /** */
  SubsetCellFilter(const CellFilter& superset,
    const CellPredicate& predicate);

  /** */
  virtual ~SubsetCellFilter();

  /** Return the dimension of the cells that will be identified
   * by this filter when acting on the given mesh */
  virtual int dimension(const Mesh& mesh) const 
    {return superset_.dimension(mesh);}

  /** */
  virtual XMLObject toXML() const ;

  /** */
  virtual bool lessThan(const CellFilterStub* other) const ;

  /** */
  virtual std::string description() const 
    {return "Subset(predicate=" + predicate_.description() 
        + ", super=" + superset_.description()+")";}


  /* */
  GET_RCP(CellFilterStub);
    

protected:
  /** */
  virtual CellSet internalGetCells(const Mesh& mesh) const ;

private:
  /** */
  CellFilter superset_;

  /** */
  CellPredicate predicate_;
};
}

#endif
