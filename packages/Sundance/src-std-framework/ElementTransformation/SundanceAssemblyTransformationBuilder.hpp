/*
 * SundanceAssemblyTransformationBuilder.hpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#ifndef SUNDANCEASSEMBLYTRANSFORMATIONBUILDER_HPP_
#define SUNDANCEASSEMBLYTRANSFORMATIONBUILDER_HPP_

#include "SundanceIntegralGroup.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceCellType.hpp"
#include "SundanceTransformationBase.hpp"
#include "SundanceTransformationHN.hpp"
#include "SundanceInequivalentElementTransformation.hpp"

namespace Sundance {

  class AssemblyTransformationBuilder {

  public:


    /** This is (mainly used from the Assemby Loop) <br>
	In case when the basis and the test space are different <br>
	An Integral group migh have many members, but the test and trial space combination is UNIQUE !!!
	Hanging node mesh with Hermit base might be tricky ... */
    AssemblyTransformationBuilder( const RCP<IntegralGroup>& group ,
				   const Array<RCP<DOFMapBase> >& rowMaps ,
				   const Array<RCP<DOFMapBase> >& colMaps ,
				   const Mesh& mesh);

    /** */
    virtual ~AssemblyTransformationBuilder();


    /** */
    void applyTransformsToAssembly( int groupIndex ,
				    int entryPerCell ,
				    CellType cellType ,
				    int cellDim,
				    CellType maxCellType ,
				    const CellJacobianBatch& JTrans,
				    const CellJacobianBatch& JVol,
				    const Array<int>& facetNum,
				    const RCP<Array<int> >& cellLIDs,
				    RCP<Array<double> >& A );


    /** return the preTransformation*/
    const RCP<TransformationBase>& getPreTransformation() const { return preTransformation_; }

    /** return the postTransformation*/
    const RCP<TransformationBase>& getPostTransformation() const { return postTransformation_; }

    /** verbosity level */
    int verb() const {return verb_;}

    /** set verbosity level */
    void setVerbosity(int verb) {verb_=verb;}

  private:

    /* verbosity */
    int verb_;

    /* number of columns of the matrix/vector which should be transformed */
    int nrCol_;

    /* number of rows of the matrix/vector which should be transformed */
    int nrRow_;

    /** The pre-transformation */
    mutable RCP<TransformationBase> preTransformation_;

    /** The post-transformation (might be the same as the pre-transformation, the RCP will take care of that) */
    mutable RCP<TransformationBase> postTransformation_;

    /** */
    int testFuncID_;

    /** */
    int unkFuncID_;

    /** */
    const DOFMapBase* _myRowDOFMap;

    /** */
    const DOFMapBase* _myColDOFMap;

    /** */
    mutable bool hasTransformation_;

    /** */
    bool hasPreTransformation_;
    bool hasPostTransformation_;

    /** */
    mutable bool onlyVectorTransformation_;
  };

}

#endif /* SUNDANCEASSEMBLYTRANSFORMATIONBUILDER_HPP_ */
