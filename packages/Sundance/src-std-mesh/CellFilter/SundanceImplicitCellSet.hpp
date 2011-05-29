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

#ifndef SUNDANCE_IMPLICITCELLSET_H
#define SUNDANCE_IMPLICITCELLSET_H

#include "SundanceDefs.hpp"
#include "SundanceCellSetBase.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceCellType.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * ImplicitCellSet is a cell set subtype where the set of cell LIDs
 * is never stored. Iteration is done by simply advancing the 
 * LID by one. 
 * 
 * @see CellFilter, CellSetBase, CellIterator 
 **/
class ImplicitCellSet : public CellSetBase
{
public:

  /** Construct with a mesh */
  ImplicitCellSet(const Mesh& mesh, int cellDim,
    const CellType& cellType);

  /** Returns an iterator pointing to the first element
   * in the set. */
  virtual CellIterator begin() const ;

  /** Returns a past-the-end iterator */
  virtual CellIterator end() const ;

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
  int maxLID_;
};
}


#endif
