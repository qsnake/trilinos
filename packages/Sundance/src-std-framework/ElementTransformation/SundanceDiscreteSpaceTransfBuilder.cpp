/*
 * SundanceDiscreteSpaceTransfBuilder.cpp
 *
 *  Created on: Mar 21, 2010
 *      Author: benk
 */


#include "SundanceDiscreteSpaceTransfBuilder.hpp"


using namespace Sundance;


DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder():
  verb_(0),
  basisSize_(0),
  hasTransformation_(false),
  transformation_()
{
}

DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder( const Mesh& mesh, const BasisArray& basis,
							const RCP<DOFMapBase>& map):
  verb_(0),
  basisSize_(0),
  hasTransformation_(false),
  transformation_()
{

  if (mesh.allowsHangingHodes())
    {
      // in this case we have to create a transformation
      const HNDoFMapBase* myMap
	= dynamic_cast<const HNDoFMapBase*>(map.get());
      if (myMap != 0)
	{
	  // create one hanging node transformation
	  SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder basis.size():" << basis.size());
	  basisSize_ = basis.size();
	  transformation_ = rcp((TransformationBase*)(new TransformationHN( myMap , 1 , basis.size() ))); //todo: basis.size is wrong
	  // the transformation is always needed when we have a mesh with hanging nodes
	  hasTransformation_ = true;
	}
    }
  else
    {
      const MixedDOFMap * myMap = dynamic_cast<const MixedDOFMap *>( map.get() );
      if (myMap != 0)
	{
	  SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder basis.size():" << basis.size());
	  basisSize_ = basis.size();
	  transformation_ = rcp((TransformationBase*)(new InequivalentElementTransformation( mesh , myMap )));
	  hasTransformation_ = transformation_->doesAnyTransformation();
	}
    }
}

void DiscreteSpaceTransfBuilder::getDoFsWithTransformation( const Array<int>& dofs,
							    const Array<int>& functionIDs ,
							    const int chunkID ,
							    const int nrDoFsPerCell ,
							    const int nrFunctions ,
							    const int cellDim ,
							    const Array<int>& cellLIDs,
							    const RCP<GhostView<double> >& ghostView ,
							    Array<double>& localValues) const 
{
  if (hasTransformation_) 
    {
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation()" );
      // get all the DoFs for this functionID
      int DoFPerElement = nrDoFsPerCell;
      Array<double> tmpArray(DoFPerElement*cellLIDs.size());
      Array<int> facetIndex(1); //just needed for the interface
      int cellI , elemDof;
      
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation()  nrFunctions:" << nrFunctions
		     << " DoFPerElement:" << DoFPerElement <<  "  localValues.size():" << localValues.size());
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation() nrDoFsPerCell:" << nrDoFsPerCell);
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation() localValues:" << localValues);
      
      for (int nf = 0 ; nf < nrFunctions ; nf++)
	{
	  // copy the element values into the temporary array
	  for(cellI = 0 ; cellI < cellLIDs.size() ; cellI++) {
	    for (elemDof = 0 ; elemDof < DoFPerElement ; elemDof++)
	      tmpArray[cellI*DoFPerElement + elemDof] = localValues[(nrFunctions*cellI+nf)*DoFPerElement + elemDof];
	  }
	   SUNDANCE_MSG2( verb() ,"getDoFsWithTransformation() before Transformation:" << tmpArray);
	   SUNDANCE_MSG2( verb() ,"getDoFsWithTransformation() chunk:" << chunkID <<" functionID:" << functionIDs[nf]);
      // make the transformation for all elements once
	   transformation_->preapplyTranspose(
			 cellDim ,
			 functionIDs[nf] ,
			 cellLIDs ,
			 facetIndex,
			 tmpArray );
	   // copy the element values back
	   SUNDANCE_MSG2( verb() , "getDoFsWithTransformation() after Transformation:" << tmpArray );

	  for(cellI = 0 ; cellI < cellLIDs.size() ; cellI++) {
	    for (elemDof = 0 ; elemDof < DoFPerElement ; elemDof++)
	      localValues[(nrFunctions*cellI+nf)*DoFPerElement + elemDof] = tmpArray[cellI*DoFPerElement + elemDof];
	  }
	}
    }
}
