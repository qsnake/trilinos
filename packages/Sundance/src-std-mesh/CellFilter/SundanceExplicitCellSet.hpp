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

#ifndef SUNDANCE_EXPLICITCELLSET_H
#define SUNDANCE_EXPLICITCELLSET_H

#include "SundanceDefs.hpp"
#include "SundanceCellSetBase.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * ExplicitCellSet is a cell set subtype where the cell LIDs
 * are stored explicitly in an STL set. 
 * 
 * @see CellFilter, CellSet, CellSetBase, CellIterator 
 **/
class ExplicitCellSet : public CellSetBase
{
public:

  /** Construct with a mesh, initializing to an empty set */
  ExplicitCellSet(const Mesh& mesh, int cellDim,
    const CellType& cellType);

  /** Construct with a set of cells */
  ExplicitCellSet(const Mesh& mesh, int cellDim,
    const CellType& cellType,
    const Set<int>& cellLIDs);

  /** Returns an iterator pointing to the first element
   * in the set. */
  virtual CellIterator begin() const ;

  /** Returns a past-the-end iterator */
  virtual CellIterator end() const ;

  /** Returns a modifiable reference to the set of cells */
  Set<int>& cells() {return cells_;}

  /** */
  bool internalLessThan(const CellSetBase* other) const ;

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const ;
  //@}

  /* Handleable interface */
  GET_RCP(CellSetBase);

private:

  /** The set of cell LIDs */
  Set<int> cells_;

      
};
}


#endif
