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

#include "SundanceDiscreteSpace.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceFunctionSupportResolver.hpp"
#include "SundanceOut.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "TSFProductVectorSpaceImpl.hpp"
#include "Teuchos_MPIContainerComm.hpp"

using namespace Sundance;
using namespace Teuchos;

const int* vecPtr(const Array<int>& x)
{
  static const int* dum = 0;
  if (x.size()==0) return dum;
  else return &(x[0]);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb), 
    map_(),
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(1), BasisArray(tuple(basis)));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(basis.size()), basis);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const Array<CellFilter>& funcDomains,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(funcDomains, basis);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
  const CellFilter& funcDomains,
  const VectorType<double>& vecType,
  int setupVerb)
  : map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(tuple(funcDomains), BasisArray(tuple(basis)));
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const CellFilter& funcDomains,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(Array<CellFilter>(basis.size(), funcDomains), basis);
}



DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const RCP<DOFMapBase>& map,
  const RCP<Array<int> >& bcIndices,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(map), 
    mesh_(mesh),
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(map->funcDomains(), basis, bcIndices, true);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const RCP<DOFMapBase>& map,
  const VectorType<double>& vecType,
  int setupVerb)
  : map_(map), 
    mesh_(mesh),
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(map->funcDomains(), basis);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
  const SpectralBasis& spBasis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(),
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(spBasis.nterms()), 
    replicate(basis, spBasis.nterms()));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const SpectralBasis& spBasis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(basis.size() * spBasis.nterms()), 
    replicate(basis, spBasis.nterms()));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const RCP<FunctionSupportResolver>& fsr,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(basis),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(new DiscreteSpaceTransfBuilder())
{
  bool partitionBCs = false;
  DOFMapBuilder builder(mesh, fsr, partitionBCs, setupVerb);

  map_ = builder.colMap()[0];
  Array<Set<CellFilter> > cf = builder.unkCellFilters()[0];

  for (int i=0; i<cf.size(); i++)
  {
    Array<Array<CellFilter> > dimCF(mesh.spatialDim()+1);
    for (Set<CellFilter>::const_iterator 
           iter=cf[i].begin(); iter != cf[i].end(); iter++)
    {
      CellFilter f = *iter;
      int dim = f.dimension(mesh);
      dimCF[dim].append(f);
    }
    for (int d=mesh.spatialDim(); d>=0; d--)
    {
      if (dimCF[d].size() == 0) continue;
      CellFilter f = dimCF[d][0];
      for (int j=1; j<dimCF[d].size(); j++)
      {
        f = f + dimCF[d][j];
      }
      subdomains_.append(f);
      break;
    }
  }
  RCP<Array<int> > dummyBCIndices;
  
  // set up the transformation
  transformationBuilder_ = rcp(new DiscreteSpaceTransfBuilder( mesh , basis , map_ ));

  initVectorSpace(dummyBCIndices, partitionBCs);
  initImporter();
}

void DiscreteSpace::init(
  const Array<CellFilter>& regions,
  const BasisArray& basis)
{
  RCP<Array<int> > dummyBCIndices;
  init(regions, basis, dummyBCIndices, false);
}

void DiscreteSpace::init(
  const Array<CellFilter>& regions,
  const BasisArray& basis,
  const RCP<Array<int> >& isBCIndex, 
  bool partitionBCs)
{
  basis_ = basis;
  subdomains_ = regions;
  Array<RCP<BasisDOFTopologyBase> > basisTop(basis.size());
  for (int b=0; b<basis.size(); b++)
  {
    basisTop[b] = rcp_dynamic_cast<BasisDOFTopologyBase>(basis[b].ptr());
  }

  if (map_.get()==0) 
  {
    Array<Set<CellFilter> > cf(regions.size());
    for (int i=0; i<regions.size(); i++) cf[i] = makeSet(regions[i]);
    DOFMapBuilder b(setupVerb_);
    map_ = b.makeMap(mesh_, basisTop, cf);
  }

  // set up the transformation
  transformationBuilder_ = rcp(new DiscreteSpaceTransfBuilder( mesh_ , basis , map_ ));

  initVectorSpace(isBCIndex, partitionBCs);

  if (!partitionBCs)
  {
    initImporter();
  }
}

void DiscreteSpace::initVectorSpace(
  const RCP<Array<int> >& isBCIndex, 
  bool partitionBCs)
{
  TEST_FOR_EXCEPTION(map_.get()==0, InternalError,
    "uninitialized map");

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  if (partitionBCs)
  {
    TEST_FOR_EXCEPT(isBCIndex.get() == 0);

    int nBCDofs = 0;
    for (int i=0; i<nDof; i++)
    {
      if ((*isBCIndex)[i]) nBCDofs++;
    }
    
    int nTotalBCDofs = nBCDofs;
    mesh().comm().allReduce(&nBCDofs, &nTotalBCDofs, 1, MPIComm::INT, MPIComm::SUM);
    int nTotalInteriorDofs = map_->numDOFs() - nTotalBCDofs;

    Array<int> interiorDofs(nDof - nBCDofs);
    Array<int> bcDofs(nBCDofs);
    int iBC = 0;
    int iIn = 0;
    for (int i=0; i<nDof; i++)
    {
      if ((*isBCIndex)[i]) bcDofs[iBC++] = lowDof+i;
      else interiorDofs[iIn++] = lowDof+i;
    }
    const int* bcDofPtr = vecPtr(bcDofs);
    const int* intDofPtr = vecPtr(interiorDofs);
    VectorSpace<double> bcSpace = vecType_.createSpace(nTotalBCDofs, nBCDofs,
      bcDofPtr, mesh().comm());
    VectorSpace<double> interiorSpace = vecType_.createSpace(nTotalInteriorDofs, nDof-nBCDofs,
      intDofPtr, mesh().comm());

    vecSpace_ = productSpace<double>(interiorSpace, bcSpace);
  }
  else
  {
    Array<int> dofs(nDof);
    for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;
    
    vecSpace_ = vecType_.createSpace(map_->numDOFs(),
      map_->numLocalDOFs(),
      &(dofs[0]), mesh().comm());
  }
}

void DiscreteSpace::initImporter()
{
  TEST_FOR_EXCEPTION(map_.get()==0, InternalError,
    "uninitialized map");
  TEST_FOR_EXCEPTION(vecSpace_.ptr().get()==0, InternalError,
    "uninitialized vector space");
  TEST_FOR_EXCEPTION(vecType_.ptr().get()==0, InternalError,
    "uninitialized vector type");
  
  RCP<Array<int> > ghostIndices = map_->ghostIndices();
  int nGhost = ghostIndices->size();
  int* ghosts = 0;
  if (nGhost!=0) ghosts = &((*ghostIndices)[0]);
  ghostImporter_ = vecType_.createGhostImporter(vecSpace_, nGhost, ghosts);
}


Array<CellFilter> DiscreteSpace::maximalRegions(int n) const
{
  CellFilter cf = new MaximalCellFilter();
  Array<CellFilter> rtn(n, cf);
  return rtn;
}



void DiscreteSpace::importGhosts(const Vector<double>& x,
  RCP<GhostView<double> >& ghostView) const
{
  ghostImporter_->importView(x, ghostView);
}
