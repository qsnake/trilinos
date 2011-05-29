#ifndef TSFLAPACKGENERALMATRIX_H
#define TSFLAPACKGENERALMATRIX_H

#include "SundanceDefs.hpp"
#include "TSFDenseSerialVector.hpp"
#include "Teuchos_Array.hpp"

namespace TSFExtended
{
  /**
   * Linear operator implemented as a LAPACK dense matrix.
   */

  class LAPACKGeneralMatrix
  {
  public:
    /** Construct with domain and range spaces, which should be
     * DenseSerialVectorSpace objects */
    LAPACKGeneralMatrix(int nRow, int nCol);

    /** Empty Ctor */
    LAPACKGeneralMatrix();

    /* added by ptb  */
    inline double& operator[](int i) {return data_[i];}
    inline const double& operator[](int i) const {return data_[i];}

    inline double& operator()(int i, int j)
      {return data_[i+nRows_*j];}
    inline const double& operator()(int i, int j) const
      {return data_[i+nRows_*j];}

    inline int getNumRows() const
    {
      return nRows_;
    }

    inline int getNumCols() const
    {
      return nCols_;
    }


    /** apply operator to dense serial vector.
     * in the range space */
    void apply(const DenseSerialVector& in,
                       DenseSerialVector& out) const ;



    /** apply inverse operator to a vector in the range space, returning
     * its preimage as a vector in the domain space. The solve is done
     * by factoring and backsolving. */
    void applyInverse(const DenseSerialVector& in,
                              DenseSerialVector& out) const ;

    /** apply adjoint operator to a vector in the domain space, returning
     * a vector in the range space. The default implementation throws an
     * exception */
    void applyAdjoint(const  DenseSerialVector& in,
                              DenseSerialVector& out) const ;

    /** apply inverse adjoint operator */
    void applyInverseAdjoint(const  DenseSerialVector& in,
                                     DenseSerialVector& out) const ;




    /** */
    void setElement(int i, int j, const double& aij) ;

    /** set all elements to zero */
    void zero() ;

    /** factor with a call to LAPACK's getrf() function */
    void factor() const ;

    /** write to a stream */
    void print(std::ostream& os) const ;



  protected:

    /** low-level matrix-vector multiply */
    void mvMult(bool transpose, const DenseSerialVector& in,
                DenseSerialVector& out) const ;

    /** low-level solve */
    void solve(bool transpose, const DenseSerialVector& in,
               DenseSerialVector& out) const ;


    int nRows_;

    int nCols_;

    DenseSerialVector data_;

    mutable Teuchos::Array<int> iPiv_;

    mutable bool isFactored_;

    mutable DenseSerialVector factorData_;

  };

}

#endif
