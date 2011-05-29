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

#ifndef SUNDANCE_CELLSETBASE_H
#define SUNDANCE_CELLSETBASE_H

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceCellIterator.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * CellSetBase is the base class for cell sets. There are two cell
 * set subtypes: ExplicitCellSet and ImplicitCellSet.
 *
 * @see CellFilter
 **/
class CellSetBase : public ObjectWithClassVerbosity<CellSetBase>,
                    public Sundance::Printable,
                    public Noncopyable,
                    public Sundance::Handleable<CellSetBase>
{
public:
  /** Construct, initializing to an empty set */
  CellSetBase(const Mesh& mesh, int cellDim,
    const CellType& cellType);

  /** Return an iterator pointing to the first element in the set */
  virtual CellIterator begin() const = 0 ;

  /** Return an iterator containing the past-the-end value */
  virtual CellIterator end() const = 0 ;

  /** Return the type of cells in this set */
  const CellType& cellType(const CellType& cellType) const
    {return cellType_;}

  /** The ID number of the mesh in which these cells exist */
  int meshID() const {return mesh_.id();}

  /** The dimension of the cells contained in this set */
  int dimension() const {return dim_;}
      
  /** The mesh in which these cells exist */
  const Mesh& mesh() const {return mesh_;}

  /** The type of the cells contained in this set */
  const CellType& cellType() const {return cellType_;}

  /** */
  bool lessThan(const CellSetBase* other) const ;

  /** */
  virtual bool internalLessThan(const CellSetBase* other) const = 0 ;

private:

  /** the mesh in which the set exists */
  Mesh mesh_;

  /** the type of cell in the set */
  CellType cellType_;

  /** the dimension of the cells in the set */
  int dim_;
      
};
}



#endif
