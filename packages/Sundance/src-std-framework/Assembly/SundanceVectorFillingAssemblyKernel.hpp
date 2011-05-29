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

#ifndef SUNDANCE_VECTORFILLINGASSEMBLYKERNEL_H
#define SUNDANCE_VECTORFILLINGASSEMBLYKERNEL_H

#include "SundanceDefs.hpp"
#include "SundanceAssemblyKernelBase.hpp"
#include "TSFLoadableVector.hpp"
#include "SundanceMapBundle.hpp"

namespace Sundance
{
using namespace Teuchos;


/**
 * VectorFillingAssemblyKernel provides a common implementation of
 * multivector-filling capabilities needed by matrix-vector assembly, 
 * multivector assembly, and functional gradient evaluation. 
 */
class VectorFillingAssemblyKernel : public AssemblyKernelBase
{
public:
  /**
   * Ctor takes several arguments:
   * \param dofMap is an array of DOFMap ptrs, one for each block 
   *
   * \param isBCIndex is an array of ptrs to arrays of ints (bools). The value 
   * (*isBCIndex[b])[d] indicates whether dof #d in block #b is or is not 
   * an essential BC dof. 
   *
   * \param lowestLocalIndex stores the lowest locally-owned DOF index for each 
   * block 
   * 
   * \param b multivector to be filled
   *
   * \param partitionBC whether dirichlet BCs are stored in a separate block
   *
   * \param verb verbosity level
   */
  VectorFillingAssemblyKernel(
    const Array<RCP<DOFMapBase> >& dofMap,
    const Array<RCP<Array<int> > >& isBCIndex,
    const Array<int>& lowestLocalIndex,
    Array<Vector<double> >& b,
    bool partitionBCs,
    int verb
    );

  /** */
  virtual ~VectorFillingAssemblyKernel(){;}

protected:

  /** 
   * Insert local vector values into a global vector
   */
  void insertLocalVectorBatch(bool isBCRqc,
    bool useCofacetCells,
    const Array<int>& funcID,  
    const Array<int>& funcBlock,
    const Array<int>& mvIndices,  
    const Array<double>& localValues) const ;

  /** 
   * Build fast lookup tables of DOFs for local cells.
   */
  void buildLocalDOFMaps(
    const RCP<StdFwkEvalMediator>& mediator,
    IntegrationCellSpecifier intCellSpec,
    const Array<Set<int> >& requiredFuncs)  ;

  

protected:
  const MapBundle& mapBundle() const {return mapBundle_;}

private:
  /** Multivector to be filled */
  Array<Vector<double> > b_;   
  /** Loadable view of the multivector to be filled */
  Array<Array<RCP<LoadableVector<double> > > > vec_;
  /** The map bundle collects together several DOF-mapping data structures */
  mutable MapBundle mapBundle_;

};


}



#endif
