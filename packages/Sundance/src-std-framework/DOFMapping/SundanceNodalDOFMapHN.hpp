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

#ifndef SUNDANCE_NODALDOFMAPHN_H
#define SUNDANCE_NODALDOFMAPHN_H

#include "SundanceDefs.hpp"
#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceHNDoFMapBase.hpp"
#include "SundanceMatrixStore.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * 
 */
class NodalDOFMapHN : public HNDoFMapBase
{
public:
  /** */
  NodalDOFMapHN(const Mesh& mesh, int nFuncs,
    const CellFilter& maxCellFilter, 
    int setupVerb);

      
  /** */
  virtual ~NodalDOFMapHN(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verb) const ;


	/**
	 * @param cellLID [in] the maxCell LID input
	 * @param funcID [in] the function ID
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] (we assume that the array is already pre-sized )*/
  void getTrafoMatrixForCell(
	    int cellLID,
	    int funcID,
	    int& trafoMatrixSize,
	    bool& doTransform,
	    Array<double>& transfMatrix ) const;

	/** Function to apply transformation for facets
	 * @param cellDim , the facet dimension
	 * @param cellLID , facet LID
	 * @param facetIndex , facet index in the maxCofacet
	 * @param funcID  [in] the function ID
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] (we assume that the array is already pre-sized )*/
  void getTrafoMatrixForFacet(
		  int cellDim,
		  int cellLID,
		  int facetIndex,
		  int funcID,
		  int& trafoMatrixSize,
		  bool& doTransform,
		  Array<double>& transfMatrix ) const;


  /** See subclass for docu */
  void getDOFsForHNCell(
		int cellDim,
		int cellLID,
        int funcID,
        Array<int>& dofs ,
        Array<double>& coefs ) const;

  /** */
  RCP<const MapStructure> mapStruct() const 
    {return structure_;}

  /** */
  int nFuncs() const {return nFuncs_;}

protected:

  void init();

  void computeOffsets(int localCount)  ;

  void shareRemoteDOFs(const Array<Array<int> >& remoteNodes);

  /** This is a temporary function which later could be out source to the Basis */
  void getPointLIDsForHN( int pointLID ,
		   int facetIndex ,
		   int maxCellIndex ,
		   Array<int>& glbLIDs ,
		   Array<double>& coefsArray,
		   Array<int>& nodeIndex );

  CellFilter maxCellFilter_;

  int dim_;

  int nFuncs_;

  int nElems_;

  int nNodes_;

  int nNodesPerElem_;

  int nFacets_;

  Array<int> elemDofs_;

  Array<int> nodeDofs_;

  RCP<const MapStructure> structure_;

  /** Is true if the cell has hanging node */
  Array<bool> hasCellHanging_;

  /** Is true if the node is hanging*/
  Array<bool> nodeIsHanging_;

  /** maps one Cell index to the point LIDs, where are global DoFs which  */
  Sundance::Map< int , Array<int> > cellsWithHangingDoF_globalDoFs_;

  /** maps one Cell index to the point LIDs, where are global DoFs which  */
  Sundance::Map< int , Array<int> > cells_To_NodeLIDs_;

  /** Maps the local (hanging) node (Point) LID to the global DoFs index */
  Sundance::Map< int , Array<int> >  hangingNodeLID_to_NodesLIDs_;

  /** Maps the local (hanging) node (Point) DoF to the global DoFs */
  Sundance::Map< int , Array<double> > hangindNodeLID_to_Coefs_;

  /** Maps the local (hanging) node (Point) DoF to the matrix index*/
  Sundance::Map< int , int > maxCellLIDwithHN_to_TrafoMatrix_;

  /** The object to store all the transformation matrixes */
  MatrixStore                          matrixStore_;

  /** */
  Array<int> facetLID_;

};

}


#endif
