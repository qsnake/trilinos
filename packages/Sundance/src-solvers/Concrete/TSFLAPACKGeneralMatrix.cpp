#include "TSFLAPACKGeneralMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

using namespace TSFExtended;
using namespace Teuchos;

LAPACKGeneralMatrix::LAPACKGeneralMatrix()
	: nRows_(0),
    nCols_(0),
    data_(),
    iPiv_(),
    isFactored_(false),
    factorData_()
{}

LAPACKGeneralMatrix::LAPACKGeneralMatrix(int nRows, int nCols)
	: nRows_(nRows),
    nCols_(nCols),
    data_(nRows*nCols),
    iPiv_(),
    isFactored_(false),
    factorData_()
{}


void LAPACKGeneralMatrix::apply(const DenseSerialVector& in,
                                DenseSerialVector& out) const
{
	mvMult(false, in, out);
}

void LAPACKGeneralMatrix::applyAdjoint(const DenseSerialVector& in,
                                       DenseSerialVector& out) const
{
	mvMult(true, in, out);
}

void LAPACKGeneralMatrix::applyInverse(const DenseSerialVector& in,
                                       DenseSerialVector& out) const
{
	solve(false, in, out);
}

void LAPACKGeneralMatrix::applyInverseAdjoint(const DenseSerialVector& in,
                                              DenseSerialVector& out) const
{
	solve(true, in, out);
}




void LAPACKGeneralMatrix::mvMult(bool transpose, const DenseSerialVector& in,
                                 DenseSerialVector& out) const
{
	const double* inPtr = &(in[0]);
	double* outPtr = &(out[0]);
	
	/* set the LAPACK transpose flag = "N" for no transpose, "T" for transpose */
	ETransp transFlag=NO_TRANS;
	if (transpose) transFlag=TRANS;
	
	// RAB & ADP : 7/10/2002 : We have fixed this!
  DenseSerialVector::blasObject().GEMV(transFlag, nRows_, nCols_, 1.0, 
                                        &(data_[0]),
                                        nRows_, inPtr, 1, 0.0, 
                                        outPtr, 1);
}


void LAPACKGeneralMatrix::solve(bool transpose, const DenseSerialVector& in,
                                DenseSerialVector& out) const
{
	int info = 0;

	const double* inPtr = &(in[0]);
	double* outPtr = &(out[0]);

	/* LAPACK overwrites the input vector argument. We copy the input 
	 * vector into the output vector, and then pass the output vector
	 * to the backsolve routine. */
	DenseSerialVector::blasObject().COPY(nRows_, inPtr, 1, outPtr, 1);

	/* factor if we haven't already done so */
	if (!isFactored_)
		{
			factor();
      isFactored_ = true;
		}
	double* dataPtr = const_cast<double*>(&(factorData_[0]));

	/* set the LAPACK transpose flag = "N" for no transpose, "T" for transpose */
	char transFlag='N';
	if (transpose) transFlag='T';;

	/* backsolve */
	int* pivPtr = const_cast<int*>(&(iPiv_[0]));
  LAPACK<int, double> lapackObj;
	lapackObj.GETRS(transFlag, nRows_, 1, dataPtr,
                  nRows_, pivPtr, outPtr, nRows_, &info);

  TEST_FOR_EXCEPTION(info != 0,
                     std::runtime_error,
                     "LAPACKGeneralMatrix backsolve failed with error code"
                     << info);
}








void LAPACKGeneralMatrix::setElement(int i, int j, const double& aij)
{
	isFactored_ = false;
	data_[nRows_*j + i] = aij;
}

void LAPACKGeneralMatrix::zero()
{
	isFactored_ = false;
	data_.zero();
}

void LAPACKGeneralMatrix::factor() const 
{
	int info = 0;

	iPiv_.resize(nRows_);

	int* pivPtr = const_cast<int*>(&(iPiv_[0]));

	factorData_ = data_;

	double* dataPtr = const_cast<double*>(&(factorData_[0]));

  LAPACK<int, double> lapackObj;

	lapackObj.GETRF(nRows_, nCols_, dataPtr, nRows_, pivPtr, &info);


  TEST_FOR_EXCEPTION(info != 0,
                     std::runtime_error,
                     "LAPACKGeneralMatrix factorization failed with error code"
                     << info);

	isFactored_ = true;
}

void LAPACKGeneralMatrix::print(std::ostream& os) const
{
	os << "LAPACK " << nRows_ << "-by-" << nCols_ << " matrix: " << std::endl;
	os << "[";
	for (int i=0; i<nRows_; i++)
		{
			os << "[";
			for (int j=0; j<nCols_; j++)
				{
					os << data_[i + nRows_*j];
					if (j < nCols_-1) os << ", ";
				}
			os << "]";
			if (i < nRows_-1) os << ", ";
		}
	os << "]";
}
