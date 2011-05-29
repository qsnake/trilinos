/*
 * SundanceDiscreteSpaceTransfBuilder.hpp
 *
 *  Created on: Mar 21, 2010
 *      Author: benk
 */

#ifndef SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_
#define SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceHNDoFMapBase.hpp"
#include "SundanceMixedDOFMap.hpp"
#include "SundanceTransformationHN.hpp"
#include "SundanceTransformationBase.hpp"
#include "SundanceInequivalentElementTransformation.hpp"

#include "TSFVectorType.hpp"
#include "TSFVectorDecl.hpp"

namespace Sundance {

  using namespace Teuchos;
  using namespace TSFExtended;

  /** This builds a transformation for one discrete space <br>
   *  The class also calls the transformation method of the created transformation object
   *  */
  class DiscreteSpaceTransfBuilder {
  public:

    /** Empty Ctor, in this case there will be no transformation*/
    DiscreteSpaceTransfBuilder();

    /** */
    DiscreteSpaceTransfBuilder( const Mesh& mesh, const BasisArray& basis,
				const RCP<DOFMapBase>& map );

    /** Dtor */
    virtual ~DiscreteSpaceTransfBuilder() {;}

    /** If there is a valid transformation then returns true, else false */
    const inline bool validTransformation() const { return hasTransformation_; }

	/** Function to make the local transformation for the discrete space
	 * @param dofs                 [in] the array with the DoF numbers
	 * @param functionIDs           [in] the functionIDs of this chunk
	 * @param chunkID              [in] function chunk ID
	 * @param nrDoFsPerCell        [in] nr DoF per Cell (per one chunk)
	 * @param nrFunctions          [in] nr of functions in this chunk
	 * @param ghostView            [in] the array where we have to get values from
	 * @param localValues          [out]  the transformed values */
	void getDoFsWithTransformation(const Array<int>& dofs,
			                             const Array<int>& functionIDs ,
			                             const int chunkID ,
			                             const int nrDoFsPerCell ,
			                             const int nrFunctions ,
			                 			 const int cellDim ,
			                 			 const Array<int>& cellLIDs,
			                             const RCP<GhostView<double> >& ghostView ,
			                             Array<double>& localValues) const;

  protected:

    /** for verbosity */
    int verb() const {return verb_;}

  private:

    /** verbosity */
    int verb_;

    /** Number of different function defined on the DoFMap*/
    int basisSize_;

    /** true if has a valid transformation */
    mutable bool hasTransformation_;

    /** The transformation for the Discrete space */
    mutable RCP<TransformationBase> transformation_;

  };

}

#endif /* SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_ */
