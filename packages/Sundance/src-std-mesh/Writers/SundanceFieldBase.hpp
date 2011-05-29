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

#ifndef SUNDANCE_FIELDBASE_H
#define SUNDANCE_FIELDBASE_H

#include "SundanceDefs.hpp"
#include "SundanceHandleable.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellSet.hpp"



namespace Sundance
{
class CellFilter;
}


namespace Sundance
{
class Mesh;
using Sundance::CellFilter;

using Sundance::Map;
using Teuchos::Array;
using Teuchos::RefCountPtr;
/**
 *
 */
class FieldBase : public Sundance::Handleable<FieldBase>
{
public:
  /** */
  FieldBase(){;}

  /** virtual dtor */
  virtual ~FieldBase(){;}

  /** */
  virtual int numElems() const {return 1;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const = 0 ;

  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const = 0 ;

  /** */
  virtual bool isPointData() const = 0 ;

  /** */
  virtual bool isCellData() const {return !isPointData();}

  /** Get a batch of data. 
   * \param batch Output array of data values. This is a 2D array packed
   * into a 1D vector with function index as the faster running index.
   */
  virtual void getDataBatch(int cellDim, const Array<int>& cellID,
    const Array<int>& funcElem, Array<double>& batch) const ;

  /**
   * Return the cell filter on which this field is defined 
   */
  virtual const CellFilter& domain() const ;

};

using Sundance::CellSet;
using Sundance::CellFilter;

CellSet connectedNodeSet(const CellFilter& f, const Mesh& mesh);
RCP<Array<int> > cellSetToLIDArray(const CellSet& cs);
}

#endif
