/*
 * SundanceInequivalentTransformation.hpp
 *
 *  Created on: 5 April, 2010
 *      Author: R.C. Kirby
 */

#ifndef SUNDANCEINEQUIVALENTELEMENTTRANSFORMATION_HPP_
#define SUNDANCEINEQUIVALENTELEMENTTRANSFORMATION_HPP_

#include "SundanceIntegralGroup.hpp"
#include "SundanceTransformationBase.hpp"
#include "SundanceMixedDOFMap.hpp"

namespace Sundance {

  class InequivalentElementTransformation : public TransformationBase {
  public:
    
    InequivalentElementTransformation( const Mesh &mesh, 
				       const MixedDOFMap *map );
    
    virtual ~InequivalentElementTransformation() {;}
    
    /** The transformation method */
    // this will potentially used in assembly  process
    virtual void preApply( const int funcID,
    		   int cellDim ,
			   const CellJacobianBatch& JTrans,
			   const CellJacobianBatch& JVol,
			   const Array<int>& facetIndex,
			   const RCP<Array<int> >& cellLIDs,
			   RCP<Array<double> >& A
			   ) const;
    
    /** */
    // this will potentially used in assembly  process
    virtual void postApply( const int funcID,
    		    int cellDim ,
			    const CellJacobianBatch& JTrans,
			    const CellJacobianBatch& JVol,
			    const Array<int>& facetIndex,
			    const RCP<Array<int> >& cellLIDs,
			    RCP<Array<double> >& A
			    ) const;

    /** */
    // this will potentially used in scatter process
    virtual void preapplyTranspose( const int cellDim,
				    const int funcID,
				    const Array<int>& cellLIDs,
				    const Array<int>& facetIndex,
				    Array<double>& A ) const;



  protected:

    /** */
    int verb() const { return verb_; }

    /** */
    void setverb(int c) { verb_ = c; }

  private :
    /** reference to the mesh, needed for preapply transpose interface */
    const Mesh &mesh_;

    /** reference to the basis; transformation calls are forwarded to it */
    const MixedDOFMap *map_;

    /** pointers to bases for each chunk */
    Array<const BasisFamilyBase *> chunkBases_;

    /** verbosity atribute */
    int verb_;
  };

}

#endif /* SUNDANCEINEQUIVALENTELEMENTTRANSFORMATION_HPP_ */
