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

#ifndef SUNDANCE_DIMENSIONALCELLFILTER_H
#define SUNDANCE_DIMENSIONALCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilterBase.hpp"


namespace Sundance
{
using namespace Teuchos;
  


/**
 * DimensionalCellFilter is a filter that identifies all cells of a 
 * specified dimension. 
 *
 * <h4> Example: </h4> get all faces in a 3D mesh
 *
 * \code
 * Mesh myMesh = myReader.getMesh();
 * CellFilter faceFilter = new DimensionalCellFilter(2);
 * CellSet faces = faceFilter.getCells(myMesh);
 * \endcode
 */
class DimensionalCellFilter : public CellFilterBase 
{
public:
  /** */
  DimensionalCellFilter(int dim);

  /** */
  virtual ~DimensionalCellFilter(){;}

  /** */
  virtual int dimension(const Mesh& mesh) const {return dim_;}

  /** */
  virtual XMLObject toXML() const ;

  /** */
  virtual std::string typeName() const {return "DimensionalCellFilter";}

  /** */
  virtual std::string description() const 
    {return "Cells(d=" + Teuchos::toString(dim_) + ")";}

  /** */
  virtual bool lessThan(const CellFilterStub* other) const ;

  /* */
  GET_RCP(CellFilterStub);

protected:
  /** get the cells */
  virtual CellSet internalGetCells(const Mesh& mesh) const ;
    
private:
    
  int dim_;


};

}



#endif
