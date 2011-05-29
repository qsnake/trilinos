/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
 /* @HEADER@ */

#ifndef TSFVECTORDECL_HPP
#define TSFVECTORDECL_HPP

#include "SundanceDefs.hpp"
#include "SundanceHandle.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFAccessibleVector.hpp"
#include "TSFRawDataAccessibleVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifdef TRILINOS_6
#include "Thyra_ProductVector.hpp"
#else
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#endif

namespace TSFExtendedOps
{
template <class Scalar, class Node1, class Node2> class LC2;
template <class Scalar, class Node> class OpTimesLC; 
/** 
 * 
 */
enum LCSign {LCAdd = 1, LCSubtract = -1};
}

namespace TSFExtended
{
  

  /** 
   * User-level vector class. 
   *
   * <h2> Creating vectors </h2>
   *
   * Ordinarily, you will never construct a Vector directly
   * from a derived type.  Rather, the createMember() method of
   * VectorSpace is used to build a vector of the appropriate
   * type, for example,
   * \code 
   * VectorType<double> vecType = new EpetraVectorType();
   * int dimension = 100;
   * VectorSpace<double> space = vecType.createSpace(dimension);
   * Vector<double> x = space.createMember(); 
   * Vector<double> y = space.createMember(); 
   * \endcode 
   * This hides from you all the ugly
   * details of creating a particular concrete type.
   *
   * You will frequently create an empty vector to be filled in later, 
   * for example,
   * \code
   * Vector<double> y;
   * \endcode
   * Note that this vector isn't just empty, it's null. Not only does 
   * it have no values assigned, it does not have a concrete type. An 
   * call a method on a null vector will result in an error. What you 
   * <it>can</it> do with a null vector is
   * <ul>
   * <li> assign another vector to it
   * \code
   * Vector<double> x = space.createVector();
   * Vector<Scalar> y;
   * y = x.copy();
   * \endcode
   * <li> assign the result of a vector operation to it
   * \code
   * Vector<Scalar> z = a*x + b*y;
   * \endcode
   */
  template <class Scalar>
  class Vector : public Sundance::Handle<Thyra::VectorBase<Scalar> >
  {
  public:
    /** \name Constructors, Destructors, and Assignment Operators */
    //@{
    HANDLE_CTORS(Vector<Scalar>, Thyra::VectorBase<Scalar>);

    /** Construct a vector from a 2-term LC */
    template<class Node1, class Node2>
    Vector(const TSFExtendedOps::LC2<Scalar, Node1, Node2>& x);

    /** Construct a vector from an operator times a linear combination */
    template<class Node>
    Vector(const TSFExtendedOps::OpTimesLC<Scalar, Node>& x);

    /** Assign a linear combination of vectors to this vector */
    template<class Node1, class Node2>
    Vector& operator=(const TSFExtendedOps::LC2<Scalar, Node1, Node2>& x);

    /** Assign a scaled linear combination to this vector */
    template<class Node>
    Vector& operator=(const TSFExtendedOps::OpTimesLC<Scalar, Node>& x);
    //@}

    /** */
    VectorSpace<Scalar> space() const 
    {return this->ptr()->space();}

    /** Return the dimension of the vector  */
    int dim() const
    {
      return this->ptr()->space()->dim();
    }
      

    /** \name ProductVector operations */
    //@{

    /** set block  */
    void setBlock(int i, const Vector<Scalar>& v);
      
    /** get block */
    Vector<Scalar> getBlock(int i) const;

    //const Vector<Scalar> getBlock(int i) const;

    /** \name Math operations */
    //@{
    /** Multiply this vector by a constant scalar factor 
     * \code
     * this = alpha * this;
     * \endcode
     */
    Vector<Scalar>& scale(const Scalar& alpha);

    /** 
     * Add a scaled vector to this vector:
     * \code
     * this = this + alpha*x 
     * \endcode
     */
    Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x);

    /** 
     * Add a scaled vector to this vector times a constant:
     * \code
     * this = gamma*this + alpha*x 
     * \endcode
     */
    Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
                           const Scalar& gamma);
    /** 
     * Add two scaled vectors to this vector times a constant:
     * \code
     * this = alpha*x + beta*y + gamma*this
     * \endcode
     */
    Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
                           const Scalar& beta, const Vector<Scalar>& y, 
                           const Scalar& gamma);

    /** 
     * Copy the values of another vector into this vector
     * \code
     * this = x
     * \endcode
     */
    Vector<Scalar>& acceptCopyOf(const Vector<Scalar>& x);

    /** 
     * Create a new vector that is a copy of this vector 
     */
    Vector<Scalar> copy() const ;

    /** 
     * Element-by-element product (Matlab dot-star operator)
     */
    Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

    /** 
     * Write the elementwise product of \f$a\f$ and \f$b\f$ into
     * <t>this.</t>
     */
    void dotStarInto(const Vector<Scalar>& a, const Vector<Scalar>& b) const ;

    /** 
     * Element-by-element division (Matlab dot-slash operator)
     */
    Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

    /** 
     * Return element-by-element reciprocal as a new vector
     */
    Vector<Scalar> reciprocal() const ;


    /** 
     * Return element-by-element absolute value as a new vector
     */
    Vector<Scalar> abs() const ;

    /** 
     * Overwrite self with element-by-element reciprocal
     */
    Vector<Scalar>& reciprocal() ;

    /** 
     * Overwrite self with element-by-element absolute value 
     */
    Vector<Scalar>& abs() ;

    /** 
     * Set all elements to a constant value
     */
    void setToConstant(const Scalar& alpha) ;

      
    /** 
     * Take dot product with another vector
     */
    Scalar dot(const Vector<Scalar>& other) const ;

    /** 
     * Overloaded operator for dot product 
     */
    Scalar operator*(const Vector<Scalar>& other) const ;

    /**
     * Compute the 1-norm of this vector
     */
    Scalar norm1() const ;

    /**
     * Compute the 2-norm of this vector
     */
    Scalar norm2() const ;

    /**
     * Compute the weighted 2-norm of this vector
     */
    Scalar norm2(const Vector<Scalar>& weights) const ;    


    /**
     * Compute the infinity-norm of this vector
     */
    Scalar normInf() const ;

    /**
     * Set all elements to zero 
     */
    void zero();


    /** Retuen the max element */
    Scalar max() const;

    /** Return the max element and the corresponding index */
    Scalar max(int& index)const;

    /** Return the max element less than bound and the corresponding index */
    Scalar max(const Scalar& bound, int& index)const;

    /** Retuen the min element */
    Scalar min()const;

    /** Return the min element and the corresponding index */
    Scalar min(int& index)const;

    /** Return the min element greater than bound and the corresponding index */
    Scalar min(const Scalar& bound, int& index)const;



    //@}


    /** \name Element loading interface */
    //@{
    /** set a single element at the given global index */
    void setElement(OrdType globalIndex, const Scalar& value) ;

    /** add to the existing value of 
     * a single element at the given global index */
    void addToElement(OrdType globalIndex, const Scalar& value) ;

    /** set a group of elements */
    void setElements(OrdType numElems, const OrdType* globalIndices, 
                     const Scalar* values) 
    {castToLoadable()->setElements(numElems, globalIndices, values);}

    /** add to a group of elements */
    void addToElements(OrdType numElems, const OrdType* globalIndices, 
                       const Scalar* values)
    {castToLoadable()->addToElements(numElems, globalIndices, values);}

    /** Do whatever finalization steps are needed by the implementation,
        for instance, synchronizing border elements. The default implementation
        * is a no-op. */
    void finalizeAssembly() {castToLoadable()->finalizeAssembly();}
    //@}

    /** \name Element access interface */
    //@{
    /** get the element at the given global index */
    Scalar getElement(OrdType globalIndex) const ;
    //{return castToAccessible()->getElement(globalIndex);}

    /** Get a batch of elements */
    void getElements(const int* globalIndices, int numElems,
      Teuchos::Array<Scalar>& elems) const 
      {castToAccessible()->getElements(globalIndices, numElems, elems);}


    /** const bracket operator  */
    const Scalar& operator[](const SequentialIterator<Scalar>& index) const ;
    
    /** non - const bracket operator  */
    Scalar& operator[](const SequentialIterator<Scalar>& index);

    //@}

    /** \name Raw data access interface */
    //@{
    /** */
    const Scalar* dataPtr() const 
    {return castToRawDataAccessible()->dataPtr();}

    /** */
    Scalar* dataPtr() 
    {return castToRawDataAccessible()->dataPtr();}
    //@}



    /** Get a stopwtach for timing vector operations */
    static RCP<Time>& opTimer()
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("Low-level vector operations");
      return rtn;
    }

    

    Vector<Scalar> eval() const {return copy();}

    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const
    {return this->ptr().get()==vec;}

    void evalInto(Vector<Scalar>& other) const {other.acceptCopyOf(*this);}

    void addInto(Vector<Scalar>& other, TSFExtendedOps::LCSign sign) const
    {
      other.update(sign, *this);
    }

    /** Describe the vector  */
    std::string description() const
    {
      const Describable* d = 
        dynamic_cast<const Describable* >(this->ptr().get());
      if (d != 0)
        {
          return d->description();
        }
      return "Vector not describable";
    }

    void print(std::ostream& os) const ;
    

  private:

    /** Cross-cast vector pointer to an accessible vector */
    const AccessibleVector<Scalar>* castToAccessible() const ;

    /** Cross-cast vector to a loadable vector */
    LoadableVector<Scalar>* castToLoadable()  ;

    /** Cross-cast vector pointer to a raw data accessible vector */
    const RawDataAccessibleVector<Scalar>* castToRawDataAccessible() const ;

    /** Cross-cast vector to a raw data accessible vector */
    RawDataAccessibleVector<Scalar>* castToRawDataAccessible();

    /** Test for valid index */
    void boundscheck(OrdType i, int dim) const ;

    /** */
    const Scalar& localElement(const OrdType& blockIndex, const OrdType& indexInBlock) const ;

    /** */
    Scalar& localElement(const OrdType& blockIndex, const OrdType& indexInBlock) ;
      
      
  };
}

template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, const TSFExtended::Vector<Scalar>& x) 
{
  x.print(os);
  return os;
}




#endif
