/*
 * SundanceMatrixStore.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: benk
 */

#include "SundanceMatrixStore.hpp"

using namespace Sundance;
using namespace Teuchos;

MatrixStore::MatrixStore() {
	nrChunk_ = 0;
	matrixLength_.resize(0);
	matrixStore_.resize(0);
}


void MatrixStore::init(int chunkNr){
	// initialize data structure
	nrChunk_ = chunkNr;
	matrixLength_.resize(chunkNr);
	matrixStore_.resize(chunkNr);
}

int MatrixStore::addMatrix(int chunkIndex, Array<double>& M){
   double L2_Vect_norm2 , tmp;
   int foundIndex = -1;
   TEST_FOR_EXCEPTION( (chunkIndex >= nrChunk_) , RuntimeError, "MatrixStore::addMatrix" );
   int nrMatrix = matrixStore_[chunkIndex].size();

   for (int ii = 0 ; ii < nrMatrix ; ii++){
	   L2_Vect_norm2 = 0.0;
	   TEST_FOR_EXCEPTION( M.size() >  matrixStore_[chunkIndex][ii].size() , RuntimeError, "MatrixStore::addMatrix" );
	   // calculate the L2 difference between the 2 matrixes
	   for (int jj = 0 ; jj < M.size() ; jj++){
		   tmp = (M[jj] - matrixStore_[chunkIndex][ii][jj]);
		   L2_Vect_norm2 += tmp*tmp;
	   }
	   // test if the two matrixes are equal, if yes then we can return the matrix Index
       if (L2_Vect_norm2 < 1e-10){
    	   foundIndex = ii;
    	   return foundIndex;
       }
   }
   // add the matrix (because this could not be found) and return
   matrixStore_[chunkIndex].append(M);
   return nrMatrix;
}

void MatrixStore::getMatrix(int chunkIndex, int matrixIndex, Array<double>& transfMatrix) const{
	// return the corresponding matrix
	TEST_FOR_EXCEPTION( (chunkIndex >= nrChunk_) || (matrixStore_[chunkIndex].size() <= matrixIndex)
			, RuntimeError, "MatrixStore::getMatrix , invalid: " << nrChunk_ << " , " << matrixIndex);
	transfMatrix = matrixStore_[chunkIndex][matrixIndex];
}
