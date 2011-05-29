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

#ifndef SUNDANCE_HOMOGENEOUSDOFMAP_H
#define SUNDANCE_HOMOGENEOUSDOFMAP_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceDOFMapBase.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * A HomogeneousDOFMap is a DOF map for the special (and common)
 * case in which every function has the same basis and is defined
 * on every cell in the mesh. 
 */
class HomogeneousDOFMap : public DOFMapBase
{
public:
  /** */
  HomogeneousDOFMap(const Mesh& mesh, 
    const BasisFamily& basis,
    int numFuncs);
                        
  /** */
  HomogeneousDOFMap(const Mesh& mesh, 
    const BasisFamily& basis,
    const Array<CellFilter>& subregions,
    int numFuncs);

  /** */
  virtual ~HomogeneousDOFMap(){;}


     

      

  /** */
  virtual void getDOFsForCellBatch(int cellDim, const Array<int>& cellLID,
    Array<int>& dofs,  
    Array<Array<int> >& funcIDs,
    Array<int>& nNodes) const ;


  /** */
  virtual void print(std::ostream& os) const ;

private:

  /** */
  void allocate(const Mesh& mesh, 
    const BasisFamily& basis,
    int numFuncs);
      
  /** */
  void buildMaximalDofTable() const ;

  /** */
  bool hasBeenAssigned(int cellDim, int cellLID) const 
    {return dofs_[cellDim][cellLID][0] != uninitializedVal();}

  /** */
  void initMap();

  /** */
  void setDOFs(int cellDim, int cellLID, 
    int& nextDOF, bool isRemote=false);

  /** */
  void shareDOFs(int cellDim,
    const Array<Array<int> >& outgoingCellRequests);

  /** */
  void computeOffsets(int dim, int localCount);

  /** */
  const Array<int>& funcIDList() const {return funcIDOnCellSet(0);}

  static int uninitializedVal() {return -1;}

  int dim_;

  Array<Array<Array<int> > > dofs_;

  mutable Array<int> maximalDofs_;

  mutable bool haveMaximalDofs_;

  Array<Array<Array<Array<int> > > > localNodePtrs_;

  Array<int> nNodesPerCell_;

  Array<int> totalNNodesPerCell_;

  Array<Array<int> > numFacets_;

  Array<Array<int> > originalFacetOrientation_;

  bool basisIsContinuous_;

      
};
}


                  


#endif
