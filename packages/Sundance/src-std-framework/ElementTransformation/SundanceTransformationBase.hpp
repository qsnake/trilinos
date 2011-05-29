/*
 * SundanceTransformationBase.hpp
 *
 *  Created on: Mar 14, 2010
 *      Author: benk
 */

#ifndef SUNDANCETRANSFORMATIONBASE_HPP_
#define SUNDANCETRANSFORMATIONBASE_HPP_

#include "SundanceIntegralGroup.hpp"

namespace Sundance {

  class TransformationBase {
  public:

    /** Trough the IntegralGroup we should have access to all information */
    TransformationBase();

    virtual ~TransformationBase();

    /** The transformation method */
    // this will potentially used in assembly  process
    virtual void preApply( const int funcID,
		       int cellDim ,
			   const CellJacobianBatch& JTrans,
			   const CellJacobianBatch& JVol,
			   const Array<int>& facetIndex,
			   const RCP<Array<int> >& cellLIDs,
			   RCP<Array<double> >& A
			   ) const = 0;

    /** */
    // this will potentially used in assembly  process
    virtual void postApply( const int funcID,
		        int cellDim ,
			    const CellJacobianBatch& JTrans,
			    const CellJacobianBatch& JVol,
			    const Array<int>& facetIndex,
			    const RCP<Array<int> >& cellLIDs,
			    RCP<Array<double> >& A
			    ) const = 0;

    /** */
    // this will potentially used in scatter process
    virtual void preapplyTranspose( const int cellDim,
				    const int funcID,
				    const Array<int>& cellLIDs,
				    const Array<int>& facetIndex,
				    Array<double>& A
				    ) const = 0;

    /** TRUE if any of the basis needs transformation and FALSE if none */
    bool doesAnyTransformation() const {return doesAnyTransformation_;}

  protected:

    /** sets the flag */
    void setDoesAnyTransformation(bool val) { doesAnyTransformation_ = val; }

    /** */
    int verb() const { return verb_; }

    /** */
    void setverb(int c) { verb_ = c; }

  private :

    /** verbosity attribute */
    int verb_;

    /** TRUE if any transformation is needed (done in the object) and FALSE if none */
    bool doesAnyTransformation_;

  };

}

#endif /* SUNDANCETRANSFORMATIONBASE_HPP_ */
