#include "SundanceInequivalentElementTransformation.hpp"

namespace Sundance
{
  InequivalentElementTransformation::InequivalentElementTransformation( const Mesh& mesh ,
									const MixedDOFMap *map ):
    mesh_( mesh ), map_( map ), chunkBases_( map->nBasisChunks() )
  {
	bool isAnyTransfReq = false;
    // extract all the bases from the dof map
    for (int i=0;i<map->nBasisChunks();i++) 
      {
	chunkBases_[i] = dynamic_cast<const BasisFamilyBase *>(map->basis(i).get());

	if (chunkBases_[i] != 0) { isAnyTransfReq = (isAnyTransfReq || chunkBases_[i]->requiresBasisTransformation()); }
      }
    // set the flag which shows if there is any transformation needed in the array of basis
    setDoesAnyTransformation(isAnyTransfReq);
  }

  void InequivalentElementTransformation::preApply( const int funcID,
		                    int cellDim ,
						    const CellJacobianBatch& JTrans,
						    const CellJacobianBatch& JVol,
						    const Array<int>& facetIndex,
						    const RCP<Array<int> >& cellLIDs,
						    RCP<Array<double> >& A
						    ) const
  
  {
	// do integration for MaxDimCells
	if ( mesh_.spatialDim() > cellDim) return;

	if (chunkBases_[map_->chunkForFuncID( funcID )]->requiresBasisTransformation())
	{
	   CellJacobianBatch JVol1;
	   mesh_.getJacobians( mesh_.spatialDim() , *(cellLIDs.get()) , JVol1 );
	   chunkBases_[map_->chunkForFuncID( funcID )]->preApplyTransformation( mesh_.cellType( mesh_.spatialDim() ) , mesh_, *cellLIDs, JVol1 , A );
	}
  }

  void InequivalentElementTransformation::postApply( const int funcID,
		                     int cellDim ,
						     const CellJacobianBatch& JTrans,
						     const CellJacobianBatch& JVol,
						     const Array<int>& facetIndex,
						     const RCP<Array<int> >& cellLIDs,
						     RCP<Array<double> >& A
						     ) const
  {
	// do integration for MaxDimCells
	if ( mesh_.spatialDim() > cellDim) return;

	if (chunkBases_[map_->chunkForFuncID( funcID )]->requiresBasisTransformation())
	{
	   CellJacobianBatch JVol1;
	   mesh_.getJacobians( mesh_.spatialDim() , *(cellLIDs.get()) , JVol1 );
	   chunkBases_[map_->chunkForFuncID( funcID )]->postApplyTransformation( mesh_.cellType( mesh_.spatialDim() ) , mesh_ , *cellLIDs, JVol1 , A );
	}
  }
  
  void InequivalentElementTransformation::preapplyTranspose( const int cellDim,
							     const int funcID,
							     const Array<int>& cellLIDs,
							     const Array<int>& facetIndex,
							     Array<double>& A ) const
  {
	// do the transformation only when it is needed so we do not have to query the Jacobians
	if (chunkBases_[map_->chunkForFuncID( funcID )]->requiresBasisTransformation())
	{
		CellJacobianBatch JVol;
		mesh_.getJacobians( mesh_.spatialDim() , cellLIDs , JVol );
		chunkBases_[map_->chunkForFuncID( funcID )]->preApplyTransformationTranspose( mesh_.cellType( mesh_.spatialDim() ) , mesh_ ,  cellLIDs, JVol , A );
	}
  }
  

}
