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

#ifndef SUNDANCE_PARTIALELEMENTDOFMAP_H
#define SUNDANCE_PARTIALELEMENTDOFMAP_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceObjectWithVerbosity.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * PartialElementDOFMap is a DOF map specialized to the case of element-based
 * DOFs on a subset of cells in the domain. All elements must have the same
 * set of functions. 
 */
class PartialElementDOFMap : public DOFMapBase
{
public:
  /** */
  PartialElementDOFMap(const Mesh& mesh, 
    const CellFilter& subdomain,
    int nFuncs,
    int setupVerb);
      
  /** */
  virtual ~PartialElementDOFMap(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verbosity) const ;

  /** */
  RCP<const Set<int> >
  allowedFuncsOnCellBatch(int cellDim,
    const Array<int>& cellLID) const ;

  /** */
  const Array<CellFilter>& funcDomains() const {return funcDomains_;}

  /** */
  virtual void print(std::ostream& os) const ;


protected:

  void init();

  void computeOffsets(int localCount)  ;

  void shareRemoteDOFs(const Array<Array<int> >& remoteElems);

  int dim_;
  int nFuncs_;
  int nElems_;
  CellFilter subdomain_;
  Array<CellFilter> funcDomains_;
  Array<int> elemDofs_;
  RCP<MapStructure> structure_;
  RCP<const Set<int> > allFuncs_;
};

}


#endif
