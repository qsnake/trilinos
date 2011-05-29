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

#ifndef SUNDANCE_DOFMAPBUILDER_H
#define SUNDANCE_DOFMAPBUILDER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceFunctionSupportResolver.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCFMeshPair.hpp"
#include "SundanceMap.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Sundance
{


/** 
 * 
 */
class DOFMapBuilder 
{
public:
  /** */
  DOFMapBuilder(int setupVerb);
  /** */
  DOFMapBuilder(const Mesh& mesh, const RCP<FunctionSupportResolver>& fsr, 
    bool findBCCols, int setupVerb);

  /** */
  const Array<RCP<DOFMapBase> >& rowMap() const {return rowMap_;}

  /** */
  const Array<RCP<DOFMapBase> >& colMap() const {return colMap_;}

  /** */
  const Array<RCP<Array<int> > >& isBCRow() const {return isBCRow_;}

  /** */
  const Array<RCP<Array<int> > >& isBCCol() const {return isBCCol_;}


  /** */
  const Array<RCP<std::set<int> > >& remoteBCCols() const 
    {return remoteBCCols_;}

  Array<Array<RCP<BasisDOFTopologyBase> > > testBasisTopologyArray() const ;

  Array<Array<RCP<BasisDOFTopologyBase> > > unkBasisTopologyArray() const ;

  Array<Array<Set<CellFilter> > > testCellFilters() const ;

  Array<Array<Set<CellFilter> > > unkCellFilters() const ;

  const Mesh& mesh() const {return mesh_;}



  RCP<DOFMapBase> makeMap(const Mesh& mesh,
    const Array<RCP<BasisDOFTopologyBase> >& basis,
    const Array<Set<CellFilter> >& filters) ;

  bool hasOmnipresentNodalMap(const Array<RCP<BasisDOFTopologyBase> >& basis,
    const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool hasCommonDomain(const Array<Set<CellFilter> >& filters) const ;

  bool hasNodalBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const ;

  bool hasCellBasis(const Array<RCP<BasisDOFTopologyBase> >& basis) const ;

  bool allFuncsAreOmnipresent(const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool isWholeDomain(const Mesh& mesh,
    int maxFilterDim,
    const Set<CellFilter>& filters) const ;

  CellFilter getMaxCellFilter(const Array<Set<CellFilter> >& filters) const ;

  static bool& allowNodalMap() {static bool rtn=true; return rtn;}

  /** */
  void extractUnkSetsFromFSR(const FunctionSupportResolver& fsr,
    Array<Set<int> >& funcSets,
    Array<CellFilter>& regions) const ;

  /** */
  void extractVarSetsFromFSR(const FunctionSupportResolver& fsr,
    Array<Set<int> >& funcSets,
    Array<CellFilter>& regions) const ;

  /** */
  const RCP<FunctionSupportResolver>& fsr() const {return fsr_;}

  /** */
  Sundance::Map<Set<int>, Set<CellFilter> > 
  buildFuncSetToCFSetMap(const Array<Set<int> >& funcSets,
    const Array<CellFilter>& regions,
    const Mesh& mesh) const ;
        
  void getSubdomainUnkFuncMatches(const FunctionSupportResolver& fsr,
    Array<Sundance::Map<CellFilter, Set<int> > >& fmap) const ;
        
  void getSubdomainVarFuncMatches(const FunctionSupportResolver& fsr,
    Array<Sundance::Map<CellFilter, Set<int> > >& fmap) const ;

  Array<Sundance::Map<Set<int>, CellFilter> > 
  funcDomains(const Mesh& mesh,
    const Sundance::Map<CellFilter, Set<int> >& fmap,
    Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellSet> >& inputToChildrenMap) const ;

  Sundance::Map<CellFilter, Set<int> > domainToFuncSetMap(const Array<Set<CellFilter> >& filters) const ;

private:

  Set<CellFilter> reduceCellFilters(const Mesh& mesh,
    const Set<CellFilter>& inputSet) const ;

  bool hasUnks() const ;

  bool unksAreOmnipresent() const ;

  bool testsAreOmnipresent() const ;

  bool regionIsMaximal(int r) const ;

  bool isSymmetric(int block) const ;

  void markBCRows(int block) ;

  void markBCCols(int block) ;

  const MPIComm& comm() const {return mesh().comm();}

  void init(bool findBCCols);

  int verb_;

  Mesh mesh_;

  RCP<FunctionSupportResolver> fsr_;

  Array<RCP<DOFMapBase> > rowMap_;

  Array<RCP<DOFMapBase> > colMap_;

  Array<RCP<Array<int> > > isBCRow_;

  Array<RCP<Array<int> > > isBCCol_;

  Array<RCP<std::set<int> > > remoteBCCols_;

};

/** \relates DOFMapBuilder */
Array<Array<BasisFamily> > testBasisArray(const RCP<FunctionSupportResolver>& fsr) ;

/** \relates DOFMapBuilder */
Array<Array<BasisFamily> > unkBasisArray(const RCP<FunctionSupportResolver>& fsr) ;

}


#endif
