/*
 * SundanceMatrixStore.hpp
 *
 *  Created on: Jun 21, 2010
 *      Author: benk
 */

#ifndef SUNDANCEMATRIXSTORE_HPP_
#define SUNDANCEMATRIXSTORE_HPP_

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"

namespace Sundance {

using namespace Teuchos;

/** Class to store the transformation matrixes globally. <br>
 * Previously we used them one per cell (with hanging nodes) */
class MatrixStore {
public:

	/** Empty ctor*/
	MatrixStore();

	/** empty dtor */
	virtual ~MatrixStore() {;}

	/** initializes the */
	void init(int chunkNr);

	/** returns the matrix index corresponding to the chunckIndex
	 * @param chunkIndex
	 * @param M , the linearized matrix */
	int addMatrix(int chunkIndex, Array<double>& M);

	/** returns the matrix */
	void getMatrix(int chunkIndex, int matrixIndex, Array<double>& transfMatrix) const;


private:

	/** nr chunck*/
	int nrChunk_;

	/** the length of the matrix for each chunck*/
	Array<int> matrixLength_;

	/** [chunckIndex][matrixIndex] -> trafo Matrix*/
	Array< Array< Array<double> > > matrixStore_;

};

}

#endif /* SUNDANCEMATRIXSTORE_HPP_ */
