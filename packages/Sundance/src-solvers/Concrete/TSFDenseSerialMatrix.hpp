#ifndef TSF_DENSE_SERIAL_MATRIX_H
#define TSF_DENSE_SERIAL_MATRIX_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "SundancePrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "TSFSimplifiedLinearOpBaseDecl.hpp"
#include "TSFSerialVectorSpace.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFSolverState.hpp"

namespace TSFExtended
{
using namespace Teuchos;

template <class T> class LinearOperator;

/**
 * Linear operator implemented as a dense matrix.
 */

class DenseSerialMatrix : public SimplifiedLinearOpBase<double>,
                          public LoadableMatrix<double>,
                          public Sundance::Printable
{
public:
  /** Construct with domain and range spaces, which should be
   * DenseSerialVectorSpace objects */
  DenseSerialMatrix(
    const RCP<const SerialVectorSpace>& domain,
    const RCP<const SerialVectorSpace>& range);

  /** Virtual dtor */
  ~DenseSerialMatrix(){;}

  /** 
   * Apply either the operator or its transpose
   */
  virtual void applyOp(const Thyra::EOpTransp M_trans,
    const Vector<double>& in,
    Vector<double> out) const ;

  /** Insert a set of elements in a row, adding to any previously
   * existing values.  The nonzero structure of the matrix must have
   * been determined at construction time. 
   *
   * @param globalRowIndex the global index of the row to which these
   * elements belong.
   * @param nElemsToInsert the number of elements being inserted in this
   * step
   * @param globalColumnIndices array of column indices. Must 
   * be nElemsToInsert in length. 
   * @param elements array of element values. Must be nElemsToInsert in
   * length
   */
  virtual void addToRow(int globalRowIndex,
    int nElemsToInsert,
    const int* globalColumnIndices,
    const double* elementValues) ;

  /** Set all elements to zero, preserving the existing structure */
  virtual void zero() ;



  /** write to a stream */
  void print(std::ostream& os) const ;

  /** 
   * \brief Return a smart pointer for the range space 
   * for <tt>this</tt> operator.
   */
  RCP< const VectorSpaceBase<double> > range() const 
    {return range_;}

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   */
  RCP< const VectorSpaceBase<double> > domain() const 
    {return domain_;}

  /** */
  const double * const dataPtr() const {return &(data_[0]);}

  /** */
  double* dataPtr() {return &(data_[0]);}

  /** */
  int numRows() const {return nRows_;}

  /** */
  int numCols() const {return nCols_;}

  /** */
  void setRow(int row, const Array<double>& rowVals);


private:

  RCP<const SerialVectorSpace> domain_;
  RCP<const SerialVectorSpace> range_;
    
  int nRows_;
  int nCols_;
  Array<double> data_;
};


/** \relates DenseSerialMatrix */
void denseSVD(const LinearOperator<double>& A,
  LinearOperator<double>& U,  
  Vector<double>& Sigma,
  LinearOperator<double>& Vt);

/** \relates DenseSerialMatrix */
SolverState<double> denseSolve(const LinearOperator<double>& A,
  const Vector<double>& b,
  Vector<double>& x);

}

#endif
