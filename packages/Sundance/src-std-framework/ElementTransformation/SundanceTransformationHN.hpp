/*
 * SundanceTransformationHN.hpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#ifndef SUNDANCETRANSFORMATIONHN_HPP_
#define SUNDANCETRANSFORMATIONHN_HPP_

#include "SundanceTransformationBase.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceNodalDOFMapHN.hpp"
#include "SundanceHNDoFMapBase.hpp"

namespace Sundance {

class TransformationHN: public Sundance::TransformationBase {
public:

	TransformationHN(const HNDoFMapBase* dofMap ,
			const int nrCol , const int nrRaw);

	virtual ~TransformationHN();

	/** */
  virtual void preApply( const int funcID ,
	         int cellDim ,
			 const CellJacobianBatch& JTrans,
			 const CellJacobianBatch& JVol,
			 const Array<int>& facetIndex,
			 const RCP<Array<int> >& cellLIDs,
			 RCP<Array<double> >& A ) const ;

	/** */
  virtual void postApply( const int funcID ,
		      int cellDim ,
			  const CellJacobianBatch& JTrans,
			  const CellJacobianBatch& JVol,
			  const Array<int>& facetIndex,
			  const RCP<Array<int> >& cellLIDs,
			  RCP<Array<double> >& A
			  ) const ;
  
  /** */
  // this will potentially used in scatter process
  virtual void preapplyTranspose( const int cellDim,
				  const int funcID,
				  const Array<int>& cellLIDs,
				  const Array<int>& facetIndex,
				  Array<double>& A ) const;

protected:

	void multiplyFromLeftWithTransp(Array<double>& M , double* A_end , const double* A_copy) const;

	void multiplyFromRight(double* A_end , Array<double>& M , const double* A_copy) const;

	void multiplyFromLeft(Array<double>& M , double* A_end , const double* A_copy ,
			              const int nrRow, const int nrCol) const;

private:

	/** The DoF Map */
	const HNDoFMapBase* dofMap_;

	/** */
	const int nrRow_;

	/** */
	const int nrCol_;

};

}

#endif /* SUNDANCETRANSFORMATIONHN_HPP_ */
