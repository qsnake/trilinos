#ifndef ANASAZI_TSF_ADAPTER_HPP
#define ANASAZI_TSF_ADAPTER_HPP

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziConfigDefs.hpp"


#include "TSFVectorImpl.hpp"
#include "TSFVectorOpsImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "SundanceTabs.hpp"

namespace Anasazi 
{
using TSFExtended::Vector;
using Teuchos::RCP;
using Teuchos::Array;

class SimpleMV
{
public:
  SimpleMV() : data_() {}

  SimpleMV(int n) : data_(rcp(new Array<Vector<double> >(n))) {}

  SimpleMV(const Array<Vector<double> >& data) 
    : data_(rcp(new Array<Vector<double> >(data.size())))
    {
      for (int i=0; i<data.size(); i++)
      {
        (*data_)[i] = data[i].copy();
      }
    }

  SimpleMV clone() const
    {
      return SimpleMV(*data_);
    }

  Vector<double>& operator[](int i) {return (*data_)[i];}

  const Vector<double>& operator[](int i) const {return (*data_)[i];}

  int size() const {return data_->size();}

  void resize(int n)
    {
      data_->resize(n);
    }
  
private:
  RCP<Array<Vector<double> > > data_;
};



inline std::ostream& operator<<(std::ostream& os, const SimpleMV& mv)
{
  os << "MV (size=" << mv.size() << ")" << std::endl;
  for (int i=0; i<mv.size(); i++)
  {
    os << "ptr=" << mv[i].ptr().get() << std::endl;
    os << mv[i] << std::endl;
  }

  return os;
}

/** */
template<>
class MultiVecTraits<double, SimpleMV >
{
public:
  typedef SimpleMV _MV;
  typedef Teuchos::ScalarTraits<double> SCT;

  static double one() {static double rtn = SCT::one(); return rtn;}
  static double zero() {static double rtn = SCT::zero(); return rtn;}

  /** \name Creation methods */
  //@{

  /**
   */
  static RCP<_MV> Clone( const  _MV & mv, const int numvecs )
    { 
      //Out::os() << "Clone(nv == " << numvecs << ")" << endl;
      TEST_FOR_EXCEPT(mv.size() <= 0);
      TEST_FOR_EXCEPT(numvecs <= 0);

      RCP<_MV> rtn = rcp(new _MV(numvecs));
      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[0].copy();
        (*rtn)[i].setToConstant(zero());
      }
      return rtn;
    }

  /**
   *
   */
  static RCP< _MV > CloneCopy( const  _MV & mv )
    { 
      //Out::os() << "CloneCopy()" << endl;
      int numvecs = mv.size();
      TEST_FOR_EXCEPT(numvecs <= 0);

      // create the new multivector
      RCP<_MV> rtn = rcp(new _MV(numvecs));
      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[i].copy();
      }
      return rtn;
    }

  /** 
      
  */
  static RCP< _MV > CloneCopy( const  _MV & mv, const std::vector<int>& index )
    { 
      //Out::os() << "CloneCopy() indexed" << endl;
      int numvecs = index.size();
      TEST_FOR_EXCEPT(numvecs <= 0);
      TEST_FOR_EXCEPT((int) index.size() > mv.size());

//      TEST_FOR_EXCEPT(detectRepeatedIndex(index));

      // create the new multivector
      RCP<  _MV  > rtn = rcp(new _MV(numvecs));

      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[index[i]].copy();
      }
      return rtn;
    }

  /**

  */      
  static RCP< _MV > CloneViewNonConst(  _MV & mv, const std::vector<int>& index )
    {
      int numvecs = index.size();
      //Out::os() << "index.size() = " << numvecs << endl;
      //Out::os() << "input size = " << mv.size() << endl;
      TEST_FOR_EXCEPT(numvecs <= 0);
      TEST_FOR_EXCEPT((int) index.size() > mv.size());

//      TEST_FOR_EXCEPT(detectRepeatedIndex(index));

      // create the new multivector
      RCP<  _MV  > rtn = rcp(new _MV(numvecs));

      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[index[i]]; // shallow copy
      }

      return rtn;
    }

  /**
   *
   */      
  static RCP<const _MV > CloneView( const _MV & mv, const std::vector<int>& index )
    {
      //Out::os() << "CloneView()" << endl;
      int numvecs = index.size();
      //Out::os() << "index size = " << numvecs << endl;
      //Out::os() << "input size = " << mv.size() << endl;
      TEST_FOR_EXCEPT(numvecs <= 0);
      TEST_FOR_EXCEPT((int) index.size() > mv.size());

//      TEST_FOR_EXCEPT(detectRepeatedIndex(index));

      // create the new multivector
      RCP<  const _MV  > rtn = rcp(new _MV(numvecs));

      for (int i=0; i<numvecs; i++)
      {
        (*(rcp_const_cast<_MV>(rtn)))[i] = mv[index[i]]; // shallow copy
      }
      return rtn;
    }

  //@}

  /** \name Attribute methods */
  //@{

  /** Obtain the vector length of \c mv. */
  static int GetVecLength( const  _MV & mv )
    {
      TEST_FOR_EXCEPT(mv.size() <= 0);
      return mv[0].space().dim(); 
    }

  /** Obtain the number of vectors in \c mv */
  static int GetNumberVecs( const  _MV & mv )
    {
      //Out::os() << "GetNumVec(" << mv.size() << ")" << endl;
      return mv.size(); 
    }

  //@}

  /** \name Update methods */
  //@{

  /*! \brief Update \c mv with \f$ \alpha A B + \beta mv \f$.
   */
  static void MvTimesMatAddMv( const double alpha, const  _MV & A, 
    const Teuchos::SerialDenseMatrix<int,double>& B, 
    const double beta,  _MV & mv )
    {
//      Out::os() << "MvTimesMatAddMv()" << endl;
      int n = B.numCols();
//      Out::os() << "B.numCols()=" << n << endl;

      TEST_FOR_EXCEPT(mv.size() != n);

      for (int j=0; j<mv.size(); j++)
      {
        Vector<double> tmp;
        if (beta==one())
        {
          tmp = mv[j].copy();
        }
        else if (beta==zero())
        {
          tmp = mv[j].copy();
          tmp.setToConstant(zero());
        }
        else
        {
          tmp = beta * mv[j];
        }
        if (alpha != zero())
        {
          for (int i=0; i<A.size(); i++)
          {
            tmp = tmp + alpha*B(i,j)*A[i];
          }
        }
        mv[j].acceptCopyOf(tmp);
      }
    }

  /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
   */
  static void MvAddMv( const double alpha, const  _MV & A, 
    const double beta,  const  _MV & B,  _MV & mv )
    { 
      TEST_FOR_EXCEPT(A.size() != B.size());
      mv.resize(A.size());
      for (int i=0; i<A.size(); i++)
      {
        if (alpha==zero() && beta != zero()) mv[i].acceptCopyOf( beta*B[i]);
        else if (beta==zero() && alpha != zero()) mv[i].acceptCopyOf(alpha*A[i]);
        else if (alpha!=zero() && beta!=zero())
          mv[i].acceptCopyOf( (alpha*A[i]) + (beta*B[i]) ) ;
        else
        {
          mv[i].acceptCopyOf(A[i].copy());
          mv[i].setToConstant(zero());
        }
      }
    }

  /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
   */
  static void MvTransMv( const double alpha, const  _MV & A, const  _MV & mv, 
    Teuchos::SerialDenseMatrix<int,double>& B )
    { 
      // Create a multivector to hold the result (m by n)
      int m = A.size();
      int n = mv.size();
//      B.shape(m, n);
      //Out::os() << "m=" << m << ", n=" << n << endl;
      for (int i=0; i<m; i++)
      {
        for (int j=0; j<n; j++)
        {
          B(i,j) = alpha * (A[i] * mv[j]);
        }
      }
    
    }

  /**
   * Dot product
  */
  static void MvDot( const  _MV & mv, const  _MV & A, std::vector<double> &b )
    {
      //Out::os() << "MvDot()" << endl;
      TEST_FOR_EXCEPT(mv.size() != A.size());
      b.resize(A.size());
      for (int i=0; i<mv.size(); i++) 
        b[i] = mv[i] * A[i];
    }

  /** Scale each element of the vectors in \c *this with \c alpha.
   */
  static void MvScale (  _MV & mv, const double alpha )
    { 
      //Out::os() << "MvScale()" << endl;
      for (int i=0; i<mv.size(); i++) mv[i].scale(alpha);
    }
    
  /** Scale each element of the \c i-th vector in \c *this with \c alpha[i].
   */
  static void MvScale (  _MV & mv, const std::vector<double>& alpha ) 
    { 
      //Out::os() << "MvScale() vector" << endl;
      for (int i=0; i<mv.size(); i++) mv[i].scale(alpha[i]);
    }
    
  //@}

  /** \name Norm method */
  //@{

  /** Compute the 2-norm of each individual vector of \c mv. */
  static void MvNorm( const  _MV & mv, 
    std::vector<Teuchos::ScalarTraits<double>::magnitudeType> &normvec )
    { 
//      Out::os() << "MvNorm()" << endl;
      normvec.resize(mv.size());
      for (int i=0; i<mv.size(); i++) 
      {
        normvec[i] = mv[i].norm2();
        //      Out::os() << "i=" << i << " |v|=" << normvec[i] << endl;
      }
      
    }

  //@}

  /** \name Initialization methods */
  //@{

  /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
   */
  static void SetBlock( const  _MV & A, const std::vector<int>& index,  _MV & mv )
    { 
      //Out::os() << "SetBlock()" << endl;
      TEST_FOR_EXCEPT(A.size() < (int) index.size());
//      mv.resize(index.size());
//      TEST_FOR_EXCEPT(detectRepeatedIndex(index));
      for (unsigned int i=0; i<index.size(); i++)
      {
        mv[index[i]].acceptCopyOf(A[i]);
      }
    }

  /*! \brief Replace the vectors in \c mv with random vectors.
   */
  static void MvRandom(  _MV & mv )
    { 
      for (int i=0; i<mv.size(); i++) randomize(mv[i]); 
    }

  /*! \brief Replace each element of the vectors in \c mv with \c alpha.
   */
  static void MvInit(  _MV & mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { 
      //Out::os() << "MvInit()" << endl;
      for (int i=0; i<mv.size(); i++) mv[i].setToConstant(alpha); 
    }

  //@}

  /** \name Print method */
  //@{

  /** Print the \c mv multi-vector to the \c os output stream. */
  static void MvPrint( const  _MV & mv, std::ostream& os )
    { 
      os << mv << endl;
    }

  //@}

  /** */
  static bool detectRepeatedIndex(const std::vector<int>& index)
    {
      std::set<int> s;
      bool rtn = false;

      for (unsigned int i=0; i<index.size(); i++)
      {
        if (s.find(index[i]) != s.end())
        {
          //Out::os() << "detected repeated index " << index[i] << endl;
          rtn = true;
        }
        s.insert(index[i]);
      }
      
      return rtn;
    }

};        


/**

*/
template <> 
class OperatorTraits < double, SimpleMV, LinearOperator<double> >
{
public:
  typedef SimpleMV _MV;  
  /**
  */    
  static void Apply ( 
    const LinearOperator< double >& Op, 
    const  _MV & x,  
    _MV & y )
    {
      //Out::os() << "Apply()" << endl;
      y.resize(x.size());
      for (int i=0; i<x.size(); i++) 
      {
//        y[i] = Op * x[i];
        y[i].acceptCopyOf(Op * x[i]);
//        Out::os() << "i=" << i << " x=" << endl;
//        Out::os() << x[i] << endl;
//        Out::os() << "i=" << i << " y=" << endl;
//        Out::os() << y[i] << endl;
//        TEST_FOR_EXCEPT(x[i].norm2() < 1.0e-12);
      }
    }
    
};



} // end of Anasazi namespace 

#endif 
