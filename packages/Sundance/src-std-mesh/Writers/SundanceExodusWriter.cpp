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

#include "SundanceExodusWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_XMLObject.hpp"

#ifdef HAVE_SUNDANCE_EXODUS 
#include "exodusII.h"
#endif


using namespace Sundance;
using namespace Teuchos;
using namespace std;


void ExodusWriter::write() const 
{
  std::string exoFile = filename();
  std::string parFile = filename();
  if (nProc() > 1) 
  {
    exoFile = exoFile + "-" + Teuchos::toString(nProc()) + "-" + Teuchos::toString(myRank());
    parFile = parFile + "-" + Teuchos::toString(nProc()) + "-" + Teuchos::toString(myRank());
  }
  exoFile = exoFile + ".exo";
  parFile = parFile + ".pxo";

  if (nProc() > 1) writeParallelInfo(parFile);
#ifdef HAVE_SUNDANCE_EXODUS
  int ws = 8;
  int exoid = ex_create(exoFile.c_str(), EX_CLOBBER, &ws, &ws);

  TEST_FOR_EXCEPTION(exoid < 0, RuntimeError, "failure to create file "
    << filename());


  Array<CellFilter> nsFilters;
  Array<int> omniNodalFuncs;
  Array<RCP<Array<int> > > funcsForNodeset;
  Array<RCP<Array<int> > > nodesForNodeset;
  Array<int> nsID;
  Array<int> nsNodesPerSet;
  Array<int> nsNodePtr;
  RCP<Array<int> > allNodes=rcp(new Array<int>());

  Array<CellFilter> blockFilters;
  Array<int> omniElemFuncs;
  Array<RCP<Array<int> > > funcsForBlock;
  Array<RCP<Array<int> > > elemsForBlock;
  Array<int> blockID;
  Array<int> nElemsPerBlock;
  Array<int> blockElemPtr;
  RCP<Array<int> > allElems=rcp(new Array<int>());

  
  findNodeSets(nsFilters, omniNodalFuncs, funcsForNodeset,
    nodesForNodeset, nsID, nsNodesPerSet, nsNodePtr, allNodes);

  findBlocks(blockFilters, omniElemFuncs, funcsForBlock,
    elemsForBlock, blockID, nElemsPerBlock, blockElemPtr, allElems);
  
  writeMesh(exoid, nsFilters, nsID, nsNodesPerSet, nsNodePtr, allNodes );
  
  writeFields(exoid, nsFilters, omniNodalFuncs, omniElemFuncs, 
    funcsForNodeset,
    nodesForNodeset, nsID);


  ex_close(exoid);
#else
  TEST_FOR_EXCEPTION(true, RuntimeError, "Exodus not enabled");
#endif
}


void ExodusWriter::offset(Array<int>& x) const
{
  for (int i=0; i<x.size(); i++) x[i]++;
}


void ExodusWriter::writeMesh(int exoid, 
  const Array<CellFilter>& nodesetFilters,
  const Array<int>& nsID,
  const Array<int>& nNodesPerSet,
  const Array<int>& nsNodePtr,
  const RCP<Array<int> >& allNodes) const
{
#ifdef HAVE_SUNDANCE_EXODUS

  int ierr = 0;

  int dim = mesh().spatialDim();
  int nElems = mesh().numCells(dim);
  
  Array<int> ssLabels = mesh().getAllLabelsForDimension(dim-1).elements();
  int numSS = 0;

  for (int ss=0; ss<ssLabels.size(); ss++) 
  {
    if (ssLabels[ss] != 0) numSS++;
  }
  
  int numNS = nsID.size();

  
  /* initialize the output file */
  int nElemBlocks = mesh().numLabels(dim);

  ierr = ex_put_init(
    exoid, 
    filename().c_str(), 
    dim,
    mesh().numCells(0), 
    nElems,
    nElemBlocks, 
    numNS,
    numSS
    );
  TEST_FOR_EXCEPT(ierr<0);


  /* write the vertices */
  int nPts = mesh().numCells(0);

  Array<double> x(nPts);
  Array<double> y(nPts);
  Array<double> z(nPts * (dim > 2));

  for (int n=0; n<nPts; n++)
  {
    const Point& P = mesh().nodePosition(n);
    x[n] = P[0];
    y[n] = P[1];
    if (dim==3) z[n] = P[2];
  }

  if (dim==2)
  {
    ierr = ex_put_coord(exoid, &(x[0]), &(y[0]), (void*) 0);
  }
  else
  {
    ierr = ex_put_coord(exoid, &(x[0]), &(y[0]), &(z[0]));
  }

  TEST_FOR_EXCEPT(ierr < 0);

  if (dim==2)
  {
    Array<std::string> cn;
    cn.append("x");
    cn.append("y");
    Array<const char*> pp;
    getCharpp(cn, pp);
    ierr = ex_put_coord_names(exoid,(char**) &(pp[0]));
  }
  else
  {
    Array<std::string> cn;
    cn.append("x");
    cn.append("y");
    cn.append("z");
    Array<const char*> pp;
    getCharpp(cn, pp);
    ierr = ex_put_coord_names(exoid, (char**)&(pp[0]));
  }

  TEST_FOR_EXCEPT(ierr < 0);

  /* write the element blocks */
  Array<int> blockLabels = mesh().getAllLabelsForDimension(dim).elements();
  int nodesPerElem = dim+1;
  std::string eType = elemType(mesh().cellType(dim));
  for (int b=0; b<blockLabels.size(); b++)
  {
    int numBlockAttr = 0;
    Array<int> blockElemLIDs;
    Array<int> nodeLIDs;
    Array<int> orient;
    mesh().getLIDsForLabel(dim, blockLabels[b], blockElemLIDs);
    int numElemsThisBlock = blockElemLIDs.size();
    mesh().getFacetLIDs(dim, blockElemLIDs, 0, nodeLIDs, orient);
    offset(nodeLIDs);
    ierr = ex_put_elem_block(
      exoid, blockLabels[b]+1, eType.c_str(), 
      numElemsThisBlock, nodesPerElem, numBlockAttr
      );

    TEST_FOR_EXCEPT(ierr < 0);
    ierr = ex_put_elem_conn(exoid, blockLabels[b]+1, &(nodeLIDs[0]));
    TEST_FOR_EXCEPT(ierr < 0);
  }


  TEST_FOR_EXCEPT(ierr < 0);
  
  /* write the side sets */
  
  for (int ss=0; ss<ssLabels.size(); ss++)
  {
    if (ssLabels[ss]==0) continue;
    Array<int> sideLIDs;
    RCP<Array<int> > elemLIDs = rcp(new Array<int>());
    RCP<Array<int> > facets = rcp(new Array<int>());
    MaximalCofacetBatch maxCofacetBatch;

    mesh().getLIDsForLabel(dim-1, ssLabels[ss], sideLIDs);
    mesh().getMaxCofacetLIDs(sideLIDs, maxCofacetBatch);
    maxCofacetBatch.getSpecifiedCofacets(0, elemLIDs, facets);

    int numSides = sideLIDs.size();
    int numDists = 0;

    offset(sideLIDs);
    offset(*elemLIDs);
    offset(*facets);

    ierr = ex_put_side_set_param(exoid, ssLabels[ss], numSides, numDists);
    ierr = ex_put_side_set(exoid, ssLabels[ss], &((*elemLIDs)[0]), &((*facets)[0]));
  }
  
  TEST_FOR_EXCEPT(ierr < 0);

  
  if (nsID.size() > 0)
  {
    /* write the node sets */
    Array<int> nsDistPerSet(nsID.size(), 0);
    Array<int> nsDistPtr(nsID.size(), 0);
    Array<int> emptyDist(1, 0);

    offset(*allNodes);
  
    ierr = ex_put_concat_node_sets( exoid,
      (int*) &(nsID[0]),
      (int*) &(nNodesPerSet[0]),
      &(nsDistPerSet[0]),
      (int*) &(nsNodePtr[0]),
      &(nsDistPtr[0]),
      &((*allNodes)[0]),
      &(emptyDist[0]));

    TEST_FOR_EXCEPT(ierr < 0);
  }
#else
  TEST_FOR_EXCEPTION(true, RuntimeError, "Exodus not enabled");
#endif
}


void ExodusWriter::writeFields(int exoid, 
  const Array<CellFilter>& nodesetFilters,
  const Array<int>& omnipresentNodalFuncs,
  const Array<int>& omnipresentElemFuncs,
  const Array<RCP<Array<int> > >& funcsForNodeset,
  const Array<RCP<Array<int> > >& nodesForNodeset,
  const Array<int>& nsID) const 
{

#ifdef HAVE_SUNDANCE_EXODUS
  int nNodalFuncs = omnipresentNodalFuncs().size();
  int nElemFuncs = omnipresentElemFuncs().size();

  int nNodesetFuncs = pointScalarFields().size() - nNodalFuncs;

  int nNodesets = funcsForNodeset.size();



  Set<int> nsFuncSet;
  Map<int, Array<int> > funcToNSMap;

  for (int i=0; i<nNodesets; i++)
  {
    const Array<int>& f = *(funcsForNodeset[i]);
    for (int j=0; j<f.size(); j++)
    {
      nsFuncSet.put(f[j]);
      if (funcToNSMap.containsKey(f[j]))
      {
        funcToNSMap[f[j]].append(i);
      }
      else
      {
        funcToNSMap.put(f[j], tuple(i));
      }
    }
  }
  Array<int> nsFuncs = nsFuncSet.elements();
  TEST_FOR_EXCEPT(nsFuncs.size() != nNodesetFuncs);

  Map<int, int > funcIDToNSFuncIndex;
  for (int i=0; i<nNodesetFuncs; i++) funcIDToNSFuncIndex.put(nsFuncs[i],i);

  Array<Array<int> > nsFuncNodesets(nsFuncs.size());
  for (int i=0; i<nNodesetFuncs; i++)
  {
    nsFuncNodesets[i] = funcToNSMap.get(nsFuncs[i]);
  }

  Array<int> nodesetFuncTruthTable(nNodesetFuncs * nNodesets, 0);
  for (int i=0; i<nNodesetFuncs; i++)
  {
    for (int j=0; j<nsFuncNodesets[i].size(); j++)
    {
      int ns = nsFuncNodesets[i][j];
      nodesetFuncTruthTable[ns*nNodesetFuncs + i] = 1;
    }

    nsFuncNodesets[i] = funcToNSMap.get(nsFuncs[i]);
  }

  

  int ierr;

  if (nNodalFuncs > 0)
  {
    ierr = ex_put_var_param(exoid, "N", nNodalFuncs);
    TEST_FOR_EXCEPT(ierr < 0);
  }

  if (nElemFuncs > 0)
  {
    ierr = ex_put_var_param(exoid, "E", nElemFuncs);
    TEST_FOR_EXCEPT(ierr < 0);
  }


  if (nNodesets > 0)
  {
    ierr = ex_put_var_param(exoid, "M", nNodesetFuncs);
    TEST_FOR_EXCEPT(ierr < 0);
    
    ierr = ex_put_nset_var_tab(exoid, nNodesets, 
      nNodesetFuncs, &(nodesetFuncTruthTable[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    
    Array<std::string> nsFuncNames(nNodesetFuncs);
    Array<const char*> nsNameP;
    
    for (int i=0; i<nNodesetFuncs; i++)
    {
      nsFuncNames[i] = pointScalarNames()[nsFuncs[i]];
    }
    getCharpp(nsFuncNames, nsNameP);  
    
    ierr = ex_put_var_names(exoid, "M", nNodesetFuncs, (char**)&(nsNameP[0]));
    TEST_FOR_EXCEPT(ierr < 0);
  }


  Array<std::string> nodalFuncNames(nNodalFuncs);
  Array<std::string> elemFuncNames(nElemFuncs);
  Array<const char*> nNameP;
  Array<const char*> eNameP;
  
  for (int i=0; i<nNodalFuncs; i++)
  {
    nodalFuncNames[i] = pointScalarNames()[omnipresentNodalFuncs[i]];
  }
  getCharpp(nodalFuncNames, nNameP);  
  
  for (int i=0; i<nElemFuncs; i++)
  {
    elemFuncNames[i] = cellScalarNames()[omnipresentElemFuncs[i]];
  }
  getCharpp(elemFuncNames, eNameP);  
  
  if (nNodalFuncs > 0)
  {
    ierr = ex_put_var_names(exoid, "N", nNodalFuncs, (char**)&(nNameP[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    
    Array<double> funcVals;
    Array<int> nodeID(mesh().numCells(0));
    for (int i=0; i<mesh().numCells(0); i++) nodeID[i]=i;
    
    for (int i=0; i<nNodalFuncs; i++)
    {
      int f = omnipresentNodalFuncs[i];
      pointScalarFields()[f]->getDataBatch(0, nodeID, tuple(f), funcVals);
      int t = 1;
      int numNodes = funcVals.size();
      ierr = ex_put_nodal_var(exoid, t, i+1, numNodes, &(funcVals[0]));
      TEST_FOR_EXCEPT(ierr < 0);
    }
    
    for (int i=0; i<nNodesetFuncs; i++)
    {
      const Array<int>& ns = nsFuncNodesets[i];
      int fid = nsFuncs[i];
      
      for (int s=0; s<ns.size(); s++)
      {
        const Array<int>& nodes = *(nodesForNodeset[ns[s]]);
        pointScalarFields()[fid]->getDataBatch(0, nodes, tuple(fid), funcVals);
        int t = 1;
        int numNodes = funcVals.size();
        int id = nsID[ns[s]];
        ierr = ex_put_nset_var(exoid, t, i+1, id, numNodes, &(funcVals[0]));
        TEST_FOR_EXCEPT(ierr < 0);
      }
    }
  }

  if (nElemFuncs > 0)
  {
    ierr = ex_put_var_names(exoid, "E", nElemFuncs, (char**)&(eNameP[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    
    Array<double> funcVals;
    int dim = mesh().spatialDim();
    Array<int> elemID(mesh().numCells(dim));
    for (int i=0; i<mesh().numCells(dim); i++) elemID[i]=i;
    
    for (int i=0; i<nElemFuncs; i++)
    {
      int f = omnipresentElemFuncs[i];
      cellScalarFields()[f]->getDataBatch(dim, elemID, tuple(f), funcVals);
      int t = 1;
      int numElems = funcVals.size();
      ierr = ex_put_elem_var(exoid, t, i+1, 1, numElems, &(funcVals[0]));
      TEST_FOR_EXCEPT(ierr < 0);
    }
  }


#else
  TEST_FOR_EXCEPTION(true, RuntimeError, "Exodus not enabled");
#endif
  
}



std::string ExodusWriter::elemType(const CellType& type) const
{
  switch(type)
  {
    case TriangleCell:
      return "TRIANGLE";
    case TetCell:
      return "TETRA";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "cell type=" << type << " cannot be used as a "
        "maximal-dimension cell in exodus");
  }
  return "NULL"; //-Wall
}



void ExodusWriter::writeParallelInfo(const std::string& parfile) const 
{
  std::ofstream os(parfile.c_str());

  int dim = mesh().spatialDim();
  int nCells = mesh().numCells(dim);
  int nPts = mesh().numCells(0);

  os << myRank() << " " << nProc() << std::endl;

  os << nPts << std::endl;
  for (int i=0; i<nPts; i++)
  {
    os << i << " " << mesh().mapLIDToGID(0,i) 
       << " " << mesh().ownerProcID(0,i) << std::endl;
  }

  os << nCells << std::endl;
  for (int c=0; c<nCells; c++)
  {
    os << c << " " << mesh().mapLIDToGID(dim,c) 
       << " " << mesh().ownerProcID(dim,c) << std::endl;
  }

  for (int i=0; i<comments().length(); i++)
  {
    os << "# " << comments()[i] << std::endl;
  }
}





void ExodusWriter::findNodeSets(
  Array<CellFilter>& nodesetFilters,
  Array<int>& omnipresentFuncs,
  Array<RCP<Array<int> > >& funcsForNodeset,
  Array<RCP<Array<int> > >& nodesForNodeset,
  Array<int>& nsID,
  Array<int>& nNodesPerSet,
  Array<int>& nsNodePtr,
  RCP<Array<int> > allNodes
  ) const 
{
  int verb = 0;

  const Array<RCP<FieldBase> >& f = pointScalarFields();
  CellFilter maximal = new MaximalCellFilter();

  nNodesPerSet.resize(0);
  nsNodePtr.resize(0);
  nsID.resize(0);

  Map<CellFilter, RCP<Array<int> > > tmp;

  for (int i=0; i<f.size(); i++)
  {
    const CellFilter& cf = f[i]->domain(); 
    if (cf==maximal) 
    {
      SUNDANCE_MSG2(verb, "function #" << i << " is defined on all nodes");
      omnipresentFuncs.append(i);
      continue;
    }
    if (!tmp.containsKey(cf))
    {
      RCP<Array<int> > a = rcp(new Array<int>());
      tmp.put(cf, a);
    }
    SUNDANCE_MSG2(verb, "function #" << i << " is defined on CF " << cf);
    tmp[cf]->append(i);
  }

  int nodesetID=1;
  int nodePtr=0;
  nodesetFilters.resize(0);
  funcsForNodeset.resize(0);
  nodesForNodeset.resize(0);

  for (Map<CellFilter, RCP<Array<int> > >::const_iterator
         i=tmp.begin(); i!=tmp.end(); i++)
  {
    const CellFilter& cf = i->first;
    nodesetFilters.append(cf);
    funcsForNodeset.append(i->second);
    RCP<Array<int> > cells 
      = cellSetToLIDArray(connectedNodeSet(cf, mesh()));
    nodesForNodeset.append(cells);
    int nn = cells->size();
    nNodesPerSet.append(nn);
    nsID.append(nodesetID++);
    nsNodePtr.append(nodePtr);
    nodePtr += nn;
  }

  SUNDANCE_MSG2(verb, "node set IDs = " << nsID);
  SUNDANCE_MSG2(verb, "num nodes = " << nNodesPerSet);
  SUNDANCE_MSG2(verb, "node set pointers = " << nsNodePtr);


  int numNodes = nodePtr;
  allNodes->resize(numNodes);

  int k=0;
  for (int i=0; i<nsID.size(); i++)
  {
    SUNDANCE_MSG2(verb, "node set " << i << " funcs = " 
      << *funcsForNodeset[i]);
    SUNDANCE_MSG2(verb, "node set " << i 
      << " nodes = " << *nodesForNodeset[i]);
    const Array<int>& myCells = *(nodesForNodeset[i]);
    for (int c=0; c<myCells.size(); c++)
    {
      (*allNodes)[k++] = myCells[c];
    }
  }

  SUNDANCE_MSG2(verb, "all nodes = " << *allNodes);
}



void ExodusWriter::findBlocks(
  Array<CellFilter>& blockFilters,
  Array<int>& omnipresentFuncs,
  Array<RCP<Array<int> > >& funcsForBlock,
  Array<RCP<Array<int> > >& elemsForBlock,
  Array<int>& blockIDs,
  Array<int>& nElemsPerBlock,
  Array<int>& elemBlockPtr,
  RCP<Array<int> > allElems
  ) const 
{
  int verb=0;
  const Array<RCP<FieldBase> >& f = cellScalarFields();
  CellFilter maximal = new MaximalCellFilter();

  nElemsPerBlock.resize(0);
  elemBlockPtr.resize(0);
  blockIDs.resize(0);

  Map<CellFilter, RCP<Array<int> > > tmp;

  for (int i=0; i<f.size(); i++)
  {
    const CellFilter& cf = f[i]->domain(); 
    if (cf==maximal) 
    {
      SUNDANCE_MSG2(verb, "function #" << i << " is defined on all nodes");
      omnipresentFuncs.append(i);
      continue;
    }
    if (!tmp.containsKey(cf))
    {
      RCP<Array<int> > a = rcp(new Array<int>());
      tmp.put(cf, a);
    }
    SUNDANCE_MSG2(verb, "function #" << i << " is defined on CF " << cf);
    tmp[cf]->append(i);
  }

  int blockID=1;
  int blockPtr=0;
  blockFilters.resize(0);
  funcsForBlock.resize(0);
  elemsForBlock.resize(0);

  for (Map<CellFilter, RCP<Array<int> > >::const_iterator
         i=tmp.begin(); i!=tmp.end(); i++)
  {
    const CellFilter& cf = i->first;
    blockFilters.append(cf);
    funcsForBlock.append(i->second);
    RCP<Array<int> > cells 
      = cellSetToLIDArray(cf.getCells(mesh()));
    elemsForBlock.append(cells);
    int nn = cells->size();
    nElemsPerBlock.append(nn);
    blockIDs.append(blockID++);
    elemBlockPtr.append(blockPtr);
    blockPtr += nn;
  }

  SUNDANCE_MSG2(verb, "block IDs = " << blockIDs);
  SUNDANCE_MSG2(verb, "num elems = " << nElemsPerBlock);
  SUNDANCE_MSG2(verb, "block pointers = " << elemBlockPtr);


  int numElems = blockPtr;
  allElems->resize(numElems);

  int k=0;
  for (int i=0; i<blockIDs.size(); i++)
  {
    SUNDANCE_MSG2(verb, "block " << i << " funcs = " 
      << *funcsForBlock[i]);
    SUNDANCE_MSG2(verb, "block " << i 
      << " elems = " << *elemsForBlock[i]);
    const Array<int>& myCells = *(elemsForBlock[i]);
    for (int c=0; c<myCells.size(); c++)
    {
      (*allElems)[k++] = myCells[c];
    }
  }

  SUNDANCE_MSG2(verb, "all elems = " << *allElems);
}

void ExodusWriter::getCharpp(const Array<std::string>& s, Array<const char*>& p) const
{
  p.resize(s.size());
  for (int i=0; i<p.size(); i++) p[i] = s[i].c_str();
}


