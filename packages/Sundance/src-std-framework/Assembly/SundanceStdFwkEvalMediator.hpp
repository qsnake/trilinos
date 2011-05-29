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

#ifndef SUNDANCE_STDFWKEVALMEDIATOR_H
#define SUNDANCE_STDFWKEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceDiscreteFunction.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * StdFwkEvalMediator evaluates mesh-dependent functions in the 
 * standard framework. A number of subtypes are supported: 
 * QuadratureEvalMediator, which does evaluation on quadrature points,  
 * and NodalEvalMediator, which does evaluation at nodal points.  */

class StdFwkEvalMediator : public AbstractEvalMediator,
                           public Sundance::Printable
{
public:
  /** */
  StdFwkEvalMediator(const Mesh& mesh, int cellDim);

  /** */
  virtual ~StdFwkEvalMediator(){;}

  /** */
  void setCellBatch(const RCP<const Array<int> >& cellLID);

  /** */
  void setIntegrationSpec(IntegrationCellSpecifier intCellSpec);


  /** Update the cell type */
  virtual void setCellType(const CellType& cellType,
    const CellType& maxCellType,
    bool isInternalBdry) ;

  /** Return the Jacobian to be used in computing the volume of cells
      being integrated. This will not necessarily be the same as the
      Jacobian used for transformations of vectors: when integrating
      derivatives over boundaries, the volume is the volume of the facet,
      while the transformations are computed on the maximal cofacets. */
  const CellJacobianBatch& JVol() const {return *JVol_;}

  /** Return the Jacobian to be used in derivative transformations. */
  const CellJacobianBatch& JTrans() const ;

  /** When evaluating derivatives on boundaries, we evaluate basis
      functions on the maximal cofacets of the boundary cells. This function
      returns the facet index, relative to the maximal cofacet, of each
      boundary cell in the batch.  */
  const Array<int>& facetIndices() const {return *facetIndices_;}

  /** */
  const Array<int>& maxCellLIDs() const {return *maxCellLIDs_;}

  /** */
  int cellDim() const {return cellDim_;}

  /** */
  int maxCellDim() const {return mesh_.spatialDim();}

  /** */
  const CellType& cellType() const {return cellType_;}

  /** */
  const CellType& maxCellType() const {return maxCellType_;}

  /** */
  const RCP<const Array<int> >& cellLID() const {return cellLID_;}

  /** */
  const RCP<Array<int> >& cofacetCellLID() const {return maxCellLIDs_;}

  /** */
  IntegrationCellSpecifier integrationCellSpec() const {return intCellSpec_;}

  /** */
  bool cofacetCellsAreReady() const {return cofacetCellsAreReady_;}

  /** */
  bool isInternalBdry() const {return isInternalBdry_;}

  /** */
  bool forbidCofacetIntegrations() const 
    {return forbidCofacetIntegrations_;}


protected:
  const Mesh& mesh() const {return mesh_;}

  Mesh& mesh() {return mesh_;}

  bool& cacheIsValid() const {return cacheIsValid_;}

  /** */
  void setupFacetTransformations() const ;

  /** */
  Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > >& fCache() const {return fCache_;}
  /** */
  Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > >& dfCache() const {return dfCache_;}
  /** */
  Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > >& localValueCache() const {return localValueCache_;}

  /** */
  Map<const DiscreteFunctionData*, RCP<const MapStructure> >& mapStructCache() const
    {return mapStructCache_;}

  /** */
  Map<const DiscreteFunctionData*, bool>& fCacheIsValid() const {return fCacheIsValid_;}
  /** */
  Map<const DiscreteFunctionData*, bool>& dfCacheIsValid() const {return dfCacheIsValid_;}
  /** */
  Map<const DiscreteFunctionData*, bool>& localValueCacheIsValid() const {return localValueCacheIsValid_;}
      
private:
  Mesh mesh_;

  int cellDim_;

  CellType cellType_;

  CellType maxCellType_;

  bool isInternalBdry_;

  bool forbidCofacetIntegrations_;

  RCP<const Array<int> > cellLID_;

  mutable IntegrationCellSpecifier intCellSpec_;

  mutable RCP<CellJacobianBatch> JVol_;

  mutable RCP<CellJacobianBatch> JTrans_;

  mutable RCP<Array<int> > facetIndices_;

  mutable RCP<Array<int> > maxCellLIDs_;

  mutable bool cofacetCellsAreReady_;

  mutable bool cacheIsValid_;

  mutable bool jCacheIsValid_;



  /** */
  mutable Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > > fCache_;
  /** */
  mutable Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > > dfCache_; 
  /** */
  mutable Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > > localValueCache_;
  /** */
  mutable Map<const DiscreteFunctionData*, RCP<const MapStructure> > mapStructCache_;

  /** */
  mutable Map<const DiscreteFunctionData*, bool> fCacheIsValid_;
  /** */
  mutable Map<const DiscreteFunctionData*, bool> dfCacheIsValid_;
  /** */
  mutable Map<const DiscreteFunctionData*, bool> localValueCacheIsValid_;

};
}


#endif
