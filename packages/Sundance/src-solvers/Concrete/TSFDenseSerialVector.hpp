#ifndef TSFDENSESERIALVECTOR_H
#define TSFDENSESERIALVECTOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"

namespace TSFExtended
{
  

  /**\ingroup DenseSerial
   * Serial vector with math operations
   */

  class DenseSerialVector : public Teuchos::Array<double>
  {
  public:
    /** {\bf Constructors, Destructors, and Assignment Operator} */
    //@{
    /** Empty ctor */
    DenseSerialVector() : Teuchos::Array<double>() {;}
    /** Create a vector of length n */
    inline DenseSerialVector(int n) : Teuchos::Array<double>(n) {;}
    /** Create a vector of length n, and fill it with value */
    inline DenseSerialVector(int n, const double& value);
    /** Create a vector of length n, assuming responsibility for a C array */
    inline DenseSerialVector(int n, const double* cArray);
    //@}


    /** {\bf some reasonably efficient math operations} */
    //@{
    /** change sign */
    void negate();

    /** vector addition with result returned through a reference argument */
    inline void add(const DenseSerialVector& other,
                    DenseSerialVector& result) const ;
    /** self-modifying vector addition */
    void add(const DenseSerialVector& other) ;

    /** vector subtraction with result returned through a reference argument */
    inline void subtract(const DenseSerialVector& other,
                         DenseSerialVector& result) const ;
    /** self-modifying vector subtraction */
    void subtract(const DenseSerialVector& other) ;

    /** daxpy (z = a*x + y) with result returned through a reference argument */
    inline void daxpy(const DenseSerialVector& other, const double& a,
                      DenseSerialVector& result) const ;
    /** self-modifying daxpy */
    void daxpy(const DenseSerialVector& other, const double& a) ;

    /** element-by-element multiplication
     * with result returned through a reference argument */
    inline void eMult(const DenseSerialVector& other,
                      DenseSerialVector& result) const ;

    /** self-modifying element-by-element multiplication */
    void eMult(const DenseSerialVector& other) ;


    /** absolute value of each element */
    void abs()  ;

    /** return the value of the max element  */
    double max() const ;

    /** return the value of the min element  */
    double min() const ;

    /** compute the matlab style ".*" operation, i.e.,
     *      this[i] = y[i] * z[i]  */
    void dotStar(const DenseSerialVector& y,
                 const DenseSerialVector& z);


    /** compute the matlab style "./" operation, i.e.,
     *      this[i] = y[i] / z[i]  */
    void dotSlash(const DenseSerialVector& y,
                  const DenseSerialVector& z);


    /** multiplication by a scalar
     * with result returned through a reference argument */
    inline void scalarMult(const double& scalar,
                           DenseSerialVector& result) const ;
    /** self-modifying multiplication by a scalar */
    void scalarMult(const double& scalar);

    /** exponentiation by a scalar
     * with result returned through a reference argument */
    inline void scalarPow(const double& scalar,
                          DenseSerialVector& result ) const ;
    /** Self-modifying exponentiation by a scalar */
    void scalarPow(const double& scalar);

    /** dot product */
    double dot(const DenseSerialVector& other) const ;

    /** dot product with self */
    double norm2Squared() const ;

    /** 2-norm */
    double norm2() const {return sqrt(norm2Squared());}

    /** sum elements */
    inline double sumElements() const ;

    /** return maximum element value */
    double maxNorm() const ;

    /** set all elements to zero */
    void zero() {setScalar(0.0);}

    /** set all elements to the given value */
    void setScalar(const double& a);
    //@}

    /** {\bf overloaded math operators; for performance reasons,
     * avoid these in performance-critical
     * code} */
    //@{
    /** unary minus */
    inline DenseSerialVector operator-() const ;
    /** reflexive addition */
    inline DenseSerialVector& operator+=(const DenseSerialVector& other);
    /** reflexive subtraction */
    inline DenseSerialVector& operator-=(const DenseSerialVector& other);
    /** reflexive scalar mult */
    inline DenseSerialVector& operator*=(const double& scalar);
    /** reflexive scalar division */
    inline DenseSerialVector& operator/=(const double& scalar);

    /** addition */
    inline DenseSerialVector operator+(const DenseSerialVector& other) const ;
    /** subtraction */
    inline DenseSerialVector operator-(const DenseSerialVector& other) const ;
    /** dot product */
    inline double operator*(const DenseSerialVector& other) const ;
    /** scalar mult */
    inline DenseSerialVector operator*(const double& scalar) const ;
    /** scalar division */
    inline DenseSerialVector operator/(const double& scalar) const ;
    //@}

    /** write a brief description to std::string */
    std::string summary() const ;

    /** a BLAS object */
    inline static const Teuchos::BLAS<int, double>& blasObject()
    {static Teuchos::BLAS<int, double> rtn; return rtn;}

  private:
    double* x() {return &(operator[](0));}
    const double* x() const {return &(operator[](0));}
    void checkLength(const DenseSerialVector& other, 
                     const std::string& funcName) const ;

  };


  inline DenseSerialVector::DenseSerialVector(int n, const double& value)
    : Teuchos::Array<double>(n, value)
  {;}

  inline DenseSerialVector::DenseSerialVector(int n, const double* cArray)
    :  Teuchos::Array<double>(n)
  {
    memcpy(x(), (const void*) cArray, (size_t) (n*sizeof(double)));
  }

  inline void DenseSerialVector::checkLength(const DenseSerialVector& other, 
                                             const std::string& funcName) const
  {
    TEST_FOR_EXCEPTION(size() != other.size(), std::runtime_error,
                       "mismatch between operands " 
                       << summary() << " and " << other.summary()
                       << " in method " << funcName);
  }

  inline void DenseSerialVector::add(const DenseSerialVector& other,
                                     DenseSerialVector& result) const
  {
    result = *this;
    result.add(other);
  }

  inline void DenseSerialVector::subtract(const DenseSerialVector& other,
                                          DenseSerialVector& result) const
  {
    result = *this;
    result.subtract(other);
  }

  inline void DenseSerialVector::daxpy(const DenseSerialVector& other,
                                       const double& a,
                                       DenseSerialVector& result) const
  {
    result = *this;
    result.daxpy(other, a);
  }

  inline void DenseSerialVector::eMult(const DenseSerialVector& other,
                                       DenseSerialVector& result) const
  {
    result = *this;
    result.eMult(other);
  }

  inline void DenseSerialVector::scalarMult(const double& a,
                                            DenseSerialVector& result) const
  {
    result = *this;
    result.scalarMult(a);
  }


  inline DenseSerialVector DenseSerialVector::operator-() const
  {
    DenseSerialVector rtn = *this;
    rtn.negate();
    return rtn;
  }

  inline DenseSerialVector& DenseSerialVector::operator+=(const DenseSerialVector& other)
  {
    add(other);
    return *this;
  }

  inline DenseSerialVector DenseSerialVector::operator+(const DenseSerialVector& other) const
  {
    DenseSerialVector rtn = *this;
    rtn.add(other);
    return rtn;
  }

  inline DenseSerialVector& DenseSerialVector::operator-=(const DenseSerialVector& other)
  {
    subtract(other);
    return *this;
  }

  inline DenseSerialVector DenseSerialVector::operator-(const DenseSerialVector& other) const
  {
    DenseSerialVector rtn = *this;
    rtn.subtract(other);
    return rtn;
  }

  inline DenseSerialVector& DenseSerialVector::operator*=(const double& scalar)
  {
    scalarMult(scalar);
    return *this;
  }

  inline DenseSerialVector& DenseSerialVector::operator/=(const double& scalar)
  {
    TEST_FOR_EXCEPTION(scalar==0, std::runtime_error, 
                       "DenseSerialVector::operator/= divide by zero");
    scalarMult(1.0/scalar);
    return *this;
  }

  inline DenseSerialVector DenseSerialVector::operator*(const double& scalar) const
  {
    DenseSerialVector rtn = *this;
    rtn.scalarMult(scalar);
    return rtn;
  }

  inline double DenseSerialVector::operator*(const DenseSerialVector& other) const
  {
    return dot(other);
  }

  inline DenseSerialVector DenseSerialVector::operator/(const double& scalar) const
  {
    DenseSerialVector rtn = *this;
    TEST_FOR_EXCEPTION(scalar==0, std::runtime_error, 
                       "DenseSerialVector::operator/ divide by zero");
    rtn.scalarMult(1.0/scalar);
    return rtn;
  }

  inline DenseSerialVector operator*(const double& scalar,
                                     const DenseSerialVector& v)
  {
    return v*scalar;
  }


}




#endif
