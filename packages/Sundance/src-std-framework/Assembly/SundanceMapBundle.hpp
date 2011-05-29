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

#ifndef SUNDANCE_MAPBUNDLE_H
#define SUNDANCE_MAPBUNDLE_H

#include "SundanceDefs.hpp"
#include "SundanceLocalDOFMap.hpp"
#include "SundanceIntegrationCellSpecifier.hpp"

namespace Sundance
{
using namespace Teuchos;

class DOFMapBase;
class StdFwkEvalMediator;


/**
 * MapBundle collects several data structures needed for DOF mapping. For
 * each variable/equation block, it contains:
 * <ul>
 * <li> A DOFMapBase object
 * <li> An array indicating whether each local index is a BC index. BC
 * indices are skipped during fill, except for EssentialBC expressions.
 * <li> The lowest local index owned by this processor.
 * </ul>
 *
 * The main reason for this class is the need to maintain in some cases
 * two local maps: one based on a set of cells, another based on a set of cofacets. 
 * In a given calculation both might be needed. The chooseMap() method is used
 * to select the correct map for a given calculation. 
 *
 * A cofacet-based map will be needed when evaluating a boundary integral involving 
 * a non-Lagrange basis or a derivative of a Lagrange basis. The most common integrals
 * are over interior cells, for which the non-cofacet map is used regardless of basis
 * (because an interior cell has no cofacets and has all information needed for DOF mapping).
 */
class MapBundle
{
public:
  /** */
  MapBundle(
    const Array<RCP<DOFMapBase> >& dofMap,
    const Array<RCP<Array<int> > >& isBCIndex,
    const Array<int>& lowestLocalIndex,
    bool partitionBCs,
    int verb
    );

  /** 
   * Build fast lookup tables of DOFs for local cells.
   */
  void buildLocalDOFMaps(
    const RCP<StdFwkEvalMediator>& mediator,
    IntegrationCellSpecifier intCellSpec,
    const Array<Set<int> >& requiredFuncs,
    int verbosity) ;

  /** 
   *
   */
  RCP<const Array<int> > workSet(int block, bool useCofacets) const ;

  /** 
   * Return the global DOF map for the b-th block
   */
  const RCP<DOFMapBase>& dofMap(int b) const {return dofMap_[b];}

  /**
   * Return the bc indicator array for the b-th block 
   */
  const RCP<Array<int> >& isBCIndex(int b) const 
    {return isBCIndex_[b];}

  /**
   * Return the lowest index owned by the b-th block 
   */
  int lowestLocalIndex(int b) const {return lowestLocalIndex_[b];}

  /** */
  int nCells() const ;

  /** 
   * Select a local DOF map according to whether cofacet cells or ordinary
   * cells should be used.
   */
  const RCP<LocalDOFMap>& chooseMap(
    int block, bool useCofacets) const ;

  /** 
   * 
   */
  const RCP<const MapStructure>& mapStruct(
    int block,
    bool useCofacetCells) const 
    {
      const RCP<const LocalDOFMap>& choice = chooseMap(block, useCofacetCells);
      return choice->mapStruct(block);
    }

  /** 
   *
   */
  const Array<int>& localDOFs(
    int block,
    bool useCofacetCells,
    int chunk) const 
    {
      const RCP<const LocalDOFMap>& choice = chooseMap(block, useCofacetCells);
      return choice->localDOFs(block)[chunk];
    }

  /** 
   *
   */
  int nNodesInChunk(
    int block,
    bool useCofacetCells,
    int chunk) const 
    {
      const RCP<const LocalDOFMap>& choice = chooseMap(block, useCofacetCells);
      return choice->nLocalNodesPerChunk(block)[chunk];
    }

  /** 
   * Return the verbosity setting
   */
  int verb() const {return verb_;}

private:
  int verb_;
  Array<RCP<DOFMapBase> > dofMap_; 
  Array<RCP<Array<int> > > isBCIndex_;
  Array<int> lowestLocalIndex_;
  
  /** localDOFMap is the compact map built w/o cofacet DOFs. Either
   * this map or the cofacet map might be used in a given calculation */
  RCP<LocalDOFMap> localDOFMap_;
  /** localDOFMap is the compact map built with cofacet DOFs */
  RCP<LocalDOFMap> cofacetLocalDOFMap_;
};



}



#endif
