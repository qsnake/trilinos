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

#ifndef TSFVECTORIMPL_HPP
#define TSFVECTORIMPL_HPP


#include "TSFVectorDecl.hpp"
#include "Thyra_SUNDIALS_Ops.hpp"
#include "TSFIndexableVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "SundancePrintable.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorSpaceImpl.hpp"
#include "TSFSequentialIteratorImpl.hpp"
#endif




namespace TSFExtended
{
using Sundance::Out;
using Sundance::Tabs;
using Sundance::Printable;

//===========================================================================
template <class Scalar> 
void Vector<Scalar>::setBlock(int i, const Vector<Scalar>& v)
{
  Thyra::DefaultProductVector<Scalar>* pv = 
    dynamic_cast<Thyra::DefaultProductVector<Scalar>* >(this->ptr().get());
  TEST_FOR_EXCEPTION(pv == 0, std::runtime_error,
    "vector is not a product vector");
  Thyra::assign(pv->getNonconstVectorBlock(i).ptr(), *(v.ptr().ptr()));
}  



//===========================================================================
template <class Scalar> 
Vector<Scalar> Vector<Scalar>::getBlock(int i) const
{
  const Thyra::DefaultProductVector<Scalar>* pv = 
    dynamic_cast <const Thyra::DefaultProductVector<Scalar>* >(this->ptr().get());
  if (pv==0) 
  {
    TEST_FOR_EXCEPTION(i != 0, std::runtime_error,
      "Nonzero block index " << i << " into a vector that is not "
      "a product vector");
    return *this;
  }
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > b = pv->getVectorBlock(i);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > bb 
    = rcp_const_cast<Thyra::VectorBase<Scalar> >(b);
  return bb;
}

//===========================================================================

template <class Scalar> 
void Vector<Scalar>::print(std::ostream& os) const 
{

  const Printable* p = 
    dynamic_cast<const Printable* >(this->ptr().get());
  if (p != 0)
  {
    p->print(os);
    return;
  }
  const Thyra::ProductMultiVectorBase<Scalar>* pv = 
    dynamic_cast <const Thyra::ProductMultiVectorBase<Scalar>* >(this->ptr().get());
  if (pv != 0)
  {
    os << "ProductVector[" << std::endl;
    for (int i=0; i<this->space().numBlocks(); i++)
    {
      os << "block=" << i << std::endl;
      this->getBlock(i).print(os);
      os << std::endl;
    }
    os << "]" << std::endl;
    return;
  }
  else
  {
    for (SequentialIterator<Scalar> i=this->space().begin(); i!=this->space().end(); i++)
    {
      os << i.globalIndex() << '\t' << (*this)[i] << std::endl;
    }

  }
}

  


//===========================================================================
template <class Scalar> inline 
const AccessibleVector<Scalar>* Vector<Scalar>::castToAccessible() const
{
  const AccessibleVector<Scalar>* av 
    = dynamic_cast<const AccessibleVector<Scalar>*>(this->ptr().get());
  TEST_FOR_EXCEPTION(av==0, std::runtime_error,
    "Attempted to cast non-accessible vector "
    << this->description() << " to an AccessibleVector");
  return av;
}

//===========================================================================
template <class Scalar> inline 
LoadableVector<Scalar>* Vector<Scalar>::castToLoadable()
{
  LoadableVector<Scalar>* lv 
    = dynamic_cast<LoadableVector<Scalar>*>(this->ptr().get());
  TEST_FOR_EXCEPTION(lv==0, std::runtime_error,
    "Attempted to cast non-loadable vector "
    << this->description() << " to a LoadableVector");
  return lv;
}


//===========================================================================
template <class Scalar> inline 
const RawDataAccessibleVector<Scalar>* Vector<Scalar>::castToRawDataAccessible() const
{
  const RawDataAccessibleVector<Scalar>* av 
    = dynamic_cast<const RawDataAccessibleVector<Scalar>*>(this->ptr().get());
  TEST_FOR_EXCEPTION(av==0, std::runtime_error,
    "Attempted to cast non-accessible vector "
    << this->description() 
    << " to an RawDataAccessibleVector");
  return av;
}

//===========================================================================
template <class Scalar> inline 
RawDataAccessibleVector<Scalar>* Vector<Scalar>::castToRawDataAccessible() 
{
  RawDataAccessibleVector<Scalar>* av 
    = dynamic_cast<RawDataAccessibleVector<Scalar>*>(this->ptr().get());
  TEST_FOR_EXCEPTION(av==0, std::runtime_error,
    "Attempted to cast non-accessible vector "
    << this->description() 
    << " to an RawDataAccessibleVector");
  return av;
}







//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::scale(const Scalar& alpha)
{
  {
    TimeMonitor t(*opTimer());
    Thyra::Vt_S(this->ptr().ptr(), alpha);
  }
  return *this;
}




//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
  const Vector<Scalar>& x)
{
  Ptr<Thyra::VectorBase<Scalar> > p = this->ptr().ptr();
  p.assert_not_null();
  const Thyra::VectorBase<Scalar>* px = x.ptr().get();
  TEST_FOR_EXCEPT(px==0);
  {
    TimeMonitor t(*opTimer());
    Thyra::Vp_StV(p, alpha, *px);
  }

  return *this;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
{
  TimeMonitor t(*opTimer());
  if (this->ptr().get()==0) 
  {
    Vector<Scalar> me = x.space().createMember();
    this->ptr() = me.ptr();
  }
  
  Thyra::assign(this->ptr().ptr(), *(x.ptr()));

  return *this;
}

template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::copy() const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    rtn.acceptCopyOf(*this);
  }
  return rtn;
}



//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::dotStar(const Vector<Scalar>& other) const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
//    Thyra::ele_wise_prod(1.0, *(this->ptr()), *(other.ptr()), rtn.ptr().get());
    Thyra::ele_wise_prod(1.0, *(this->ptr()), *(other.ptr()), rtn.ptr().ptr());
  }
  return rtn;
}



//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::dotStarInto(
  const Vector<Scalar>& a, const Vector<Scalar>& b) const 
{
  TimeMonitor t(*opTimer());
//  Thyra::ele_wise_prod(1.0, *(a.ptr()), *(b.ptr()), this->ptr().get());
  Thyra::ele_wise_prod(1.0, *(a.ptr()), *(b.ptr()), this->ptr().ptr());
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::dotSlash(const Vector<Scalar>& other) const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    Thyra::ele_wise_divide(1.0, *(this->ptr)(), *(other.ptr()), rtn.ptr().ptr());
  }
  return rtn;
}






//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::abs() const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    rtn.acceptCopyOf(*this);
    rtn.abs();
  }
  return rtn;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::reciprocal() const 
{
  Vector<Scalar> rtn = space().createMember();
  {
    TimeMonitor t(*opTimer());
    rtn.acceptCopyOf(*this);
    rtn.reciprocal();
  }
  return rtn;
}



//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::abs()
{
  {
    TimeMonitor t(*opTimer());
    Thyra::abs(this->ptr().ptr(), *(this->ptr()));
  }
  return *this;
}
  




//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::reciprocal()
{
  TimeMonitor t(*opTimer());
  Thyra::reciprocal(this->ptr().ptr(), *(this->ptr()));

  return *this;
}

  

//===========================================================================
template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
  const Vector<Scalar>& x, 
  const Scalar& gamma)
{
  TimeMonitor t(*opTimer());

  ArrayView<const Scalar> a(&alpha, 1);
  Teuchos::Ptr<const Thyra::VectorBase<Scalar> > px = x.ptr().ptr();
  ArrayView<const Teuchos::Ptr<const Thyra::VectorBase<Scalar> > > apx(&px,1);
  
  Thyra::linear_combination(a, apx, gamma, this->ptr().ptr());

  return *this;
}




//===========================================================================
template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
  const Vector<Scalar>& x, 
  const Scalar& beta, 
  const Vector<Scalar>& y, 
  const Scalar& gamma)
{
  TimeMonitor t(*opTimer());

  Scalar a[2];
  a[0] = alpha;
  a[1] = beta;
  ArrayView<const Scalar> av(a,2);

  Ptr<const Thyra::VectorBase<Scalar> > vecs[2];
  vecs[0] = x.ptr().ptr().getConst();
  vecs[1] = y.ptr().ptr().getConst();
  ArrayView<const Ptr<const Thyra::VectorBase<Scalar> > > vv(vecs,2);

  Thyra::linear_combination(av, vv, gamma, this->ptr().ptr());

  return *this;
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::dot(const Vector<Scalar>& other) const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::dot(*(this->ptr)(), *(other.ptr()));
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::operator*(const Vector<Scalar>& other) const 
{
  return dot(other);
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm1() const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::norm_1(*(this->ptr)());
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2() const 
{
  TimeMonitor t(*opTimer());
  return Thyra::norm_2(*(this->ptr)());
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2(const Vector<Scalar>& weights) const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::norm_2(*(weights.ptr()), *(this->ptr)());
}





//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::normInf() const 
{
  TimeMonitor t(*opTimer());
    
  return Thyra::norm_inf(*(this->ptr)());
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::zero()
{
  TimeMonitor t(*opTimer());
    
  Thyra::assign(this->ptr().ptr(), 0.0);
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::setToConstant(const Scalar& alpha)
{
  TimeMonitor t(*opTimer());
    
  Thyra::assign(this->ptr().ptr(), alpha);
}



  

//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max()const
{
  TimeMonitor t(*opTimer());
  return Thyra::max(*(this->ptr)());
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max(int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar maxEl;
  Ptr<Scalar> maxElP(&maxEl);
  OrdType loc_index = -1;
  Ptr<OrdType> pind(&loc_index);
  Thyra::max(*(this->ptr()), maxElP, pind); 
  index = loc_index;
  return maxEl;
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max(const Scalar& bound, int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar maxEl;
  Ptr<Scalar> maxElP(&maxEl);
  OrdType loc_index = -1;
  Ptr<OrdType> pind(&loc_index);
  Thyra::maxLessThanBound(*(this->ptr()), bound, maxElP, pind); 
  index = loc_index;
  return maxEl;

}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min()const
{
  TimeMonitor t(*opTimer());
  return Thyra::min(*(this->ptr()));
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min(int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar minEl;
  Ptr<Scalar> minElP(&minEl);
  OrdType loc_index = -1;
  Ptr<OrdType> pind(&loc_index);
  Thyra::min(*(this->ptr()), minElP, pind); 
  index = loc_index;
  return minEl;
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min(const Scalar& bound, int& index)const
{
  TimeMonitor t(*opTimer());
  Scalar minEl;
  Ptr<Scalar> minElP(&minEl);
  OrdType loc_index = -1;
  Ptr<OrdType> pind(&loc_index);

  Thyra::minGreaterThanBound(*(this->ptr)(), bound, minElP, pind); 
  index = loc_index;
  return minEl;
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::getElement(OrdType globalIndex) const
{ 
  Thyra::ProductVectorBase<Scalar>* p 
    = dynamic_cast<Thyra::ProductVectorBase<Scalar>*>(&*this->ptr());

  if (p)
  {
    const Thyra::ProductVectorSpaceBase<Scalar>* pvs 
      = dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>*>(space().ptr().get());
    TEST_FOR_EXCEPT(pvs == 0);

    const Thyra::DefaultProductVectorSpace<Scalar>* dpvs 
      = dynamic_cast<const Thyra::DefaultProductVectorSpace<Scalar>*>(pvs);
    if (dpvs)
    {
      int blockIndex=-1;
      OrdType globOffsetInBlock=-1;
      dpvs->getVecSpcPoss(globalIndex, &blockIndex, &globOffsetInBlock);
      RCP<Thyra::VectorBase<Scalar> > vec_i 
        = p->getNonconstVectorBlock(blockIndex);
      Vector<Scalar> vv(vec_i);
      return vv.getElement(globOffsetInBlock);
    }
    else
    {
      int k = 0;
      for (int i = 0; i < pvs->numBlocks(); i++)
      {
        RCP<Thyra::VectorBase<Scalar> > vec_i 
          = p->getNonconstVectorBlock(i);
        int len = vec_i->space()->dim();
        if (globalIndex < k + len )
        {
          Vector<Scalar> vv(vec_i);
          int globalIndexWithinBlock = globalIndex - k;
          return vv.getElement(globalIndexWithinBlock);
          break;
        }
        k += len;
      }
    }
    TEST_FOR_EXCEPT(true); // should never happen
    return 0.0; // -Wall
  }
  else
  {
    Thyra::DefaultSpmdVector<Scalar>* dsv
      = dynamic_cast<Thyra::DefaultSpmdVector<Scalar>*>(this->ptr().get());
      
    if (dsv)
    {
      OrdType stride = dsv->getStride();
      OrdType low = dsv->spmdSpace()->localOffset();
      OrdType subdim = dsv->spmdSpace()->localSubDim();
      TEST_FOR_EXCEPTION( globalIndex < low || globalIndex >= low+subdim, 
        std::runtime_error,
        "Bounds violation: " << globalIndex << "is out of range [low" 
        << ", " <<  low+subdim << "]");
      return dsv->getPtr()[stride*(globalIndex - low)];
    }
    else
    {
      return castToAccessible()->getElement(globalIndex);
    }
  }
}

//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::setElement(OrdType globalIndex, const Scalar& value)
{ 
  Thyra::ProductVectorBase<Scalar>* p 
    = dynamic_cast<Thyra::ProductVectorBase<Scalar>*>(&*this->ptr());

  if (p)
  {
    const Thyra::ProductVectorSpaceBase<Scalar>* pvs 
      = dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>*>(space().ptr().get());
    TEST_FOR_EXCEPT(pvs == 0);

    const Thyra::DefaultProductVectorSpace<Scalar>* dpvs 
      = dynamic_cast<const Thyra::DefaultProductVectorSpace<Scalar>*>(pvs);
    if (dpvs)
    {
      int blockIndex=-1;
      OrdType globOffsetInBlock=-1;
      dpvs->getVecSpcPoss(globalIndex, &blockIndex, &globOffsetInBlock);
      RCP<Thyra::VectorBase<Scalar> > vec_i 
        = p->getNonconstVectorBlock(blockIndex);
      Vector<Scalar> vv(vec_i);
      vv.setElement(globOffsetInBlock, value);
    }
    else
    {
      int k = 0;
      for (int i = 0; i < pvs->numBlocks(); i++)
      {
        RCP<Thyra::VectorBase<Scalar> > vec_i 
          = p->getNonconstVectorBlock(i);
        int len = vec_i->space()->dim();
        if (globalIndex < k + len )
        {
          Vector<Scalar> vv(vec_i);
          int globalIndexWithinBlock = globalIndex - k;
          vv.setElement(globalIndexWithinBlock, value);
          break;
        }
        k += len;
      }
    }

  }
  else
  {
    Thyra::DefaultSpmdVector<Scalar>* dsv
      = dynamic_cast<Thyra::DefaultSpmdVector<Scalar>*>(this->ptr().get());
      
    if (dsv)
    {
      OrdType stride = dsv->getStride();
      OrdType low = dsv->spmdSpace()->localOffset();
      OrdType subdim = dsv->spmdSpace()->localSubDim();
      TEST_FOR_EXCEPTION( globalIndex < low || globalIndex >= low+subdim, 
        std::runtime_error,
        "Bounds violation: " << globalIndex << "is out of range [low" 
        << ", " <<  low+subdim << "]");
      dsv->getPtr()[stride*(globalIndex - low)] = value;
    }
    else
    {
      castToLoadable()->setElement(globalIndex, value);
    }
  }
}


//===========================================================================
template <class Scalar> inline 
Scalar& Vector<Scalar>::operator[](const SequentialIterator<Scalar>& iter)
{
  const OrdType& blockIndex = iter.blockIndex();
  const OrdType& indexInBlock = iter.indexInBlock();
  
  return localElement(blockIndex, indexInBlock);
} 



//===========================================================================
template <class Scalar> inline 
const Scalar& Vector<Scalar>::operator[](const SequentialIterator<Scalar>& iter) const
{
  const OrdType& blockIndex = iter.blockIndex();
  const OrdType& indexInBlock = iter.indexInBlock();
  
  return localElement(blockIndex, indexInBlock);
} 


//===========================================================================
template <class Scalar> inline 
const Scalar& Vector<Scalar>::localElement(const OrdType& blockIndex, const OrdType& indexInBlock) const
{
  const Thyra::ProductVectorBase<Scalar>* p 
    = dynamic_cast<const Thyra::ProductVectorBase<Scalar>*>(&*this->ptr());
  RCP<const Thyra::VectorBase<Scalar> > vec;
  if (p)
  {
    vec = p->getVectorBlock(blockIndex);
  }
  else
  {
    vec = this->ptr();
  }
  const Thyra::DefaultSpmdVector<Scalar>* dsv 
    = dynamic_cast<const Thyra::DefaultSpmdVector<Scalar>*>(this->ptr().get());
  if (dsv != 0)
  {
    return dsv->getPtr()[indexInBlock*dsv->getStride()];
  }
  return castToRawDataAccessible()->dataPtr()[indexInBlock];
} 



//===========================================================================
template <class Scalar> inline 
Scalar& Vector<Scalar>::localElement(const OrdType& blockIndex, const OrdType& indexInBlock) 
{
  Thyra::ProductVectorBase<Scalar>* p 
    = dynamic_cast<Thyra::ProductVectorBase<Scalar>*>(&*this->ptr());
  RCP<Thyra::VectorBase<Scalar> > vec;
  if (p)
  {
    vec = p->getNonconstVectorBlock(blockIndex);
  }
  else
  {
    vec = this->ptr();
  }
  
  Thyra::DefaultSpmdVector<Scalar>* dsv 
    = dynamic_cast<Thyra::DefaultSpmdVector<Scalar>*>(this->ptr().get());
  if (dsv != 0)
  {
    return dsv->getPtr()[indexInBlock*dsv->getStride()];
  }
  TSFExtended::RawDataAccessibleVector<Scalar>* d 
    = this->castToRawDataAccessible();
  return d->dataPtr()[indexInBlock];
} 





//#ifdef OBSOLETE_CODE


//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::addToElement(OrdType globalIndex, const Scalar& value)
{
  Thyra::ProductVectorBase<Scalar>* p 
    = dynamic_cast<Thyra::ProductVectorBase<Scalar>*>(&*this->ptr());

  if (p)
  {
    const Thyra::ProductVectorSpaceBase<Scalar>* pvs 
      = dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>*>(space().ptr().get());
    TEST_FOR_EXCEPT(pvs == 0);

    const Thyra::DefaultProductVectorSpace<Scalar>* dpvs 
      = dynamic_cast<const Thyra::DefaultProductVectorSpace<Scalar>*>(pvs);
    if (dpvs)
    {
      int blockIndex=-1;
      OrdType globOffsetInBlock=-1;
      dpvs->getVecSpcPoss(globalIndex, &blockIndex, &globOffsetInBlock);
      RCP<Thyra::VectorBase<Scalar> > vec_i 
        = p->getNonconstVectorBlock(blockIndex);
      Vector<Scalar> vv(vec_i);
      vv.addToElement(globOffsetInBlock, value);
    }
    else
    {
      int k = 0;
      for (int i = 0; i < pvs->numBlocks(); i++)
      {
        RCP<Thyra::VectorBase<Scalar> > vec_i 
          = p->getNonconstVectorBlock(i);
        int len = vec_i->space()->dim();
        if (globalIndex < k + len )
        {
          Vector<Scalar> vv(vec_i);
          int globalIndexWithinBlock = globalIndex - k;
          vv.addToElement(globalIndexWithinBlock, value);
          break;
        }
        k += len;
      }
    }

  }
  else
  {
    Thyra::DefaultSpmdVector<Scalar>* dsv
      = dynamic_cast<Thyra::DefaultSpmdVector<Scalar>*>(this->ptr().get());
      
    LoadableVector<Scalar>* loadable = dynamic_cast<LoadableVector<Scalar>*>(this->ptr().get());
    if (loadable)
    {
      loadable->addToElement(globalIndex, value);
    }
    else if (dsv)
    {
      OrdType stride = dsv->getStride();
      OrdType low = dsv->spmdSpace()->localOffset();
      OrdType subdim = dsv->spmdSpace()->localSubDim();
      TEST_FOR_EXCEPTION( globalIndex < low || globalIndex >= low+subdim, 
        std::runtime_error,
        "Bounds violation: " << globalIndex << "is out of range [low" 
        << ", " <<  low+subdim << "]");
      dsv->getPtr()[stride*(globalIndex - low)] += value;
    }
    else
    {
      TEST_FOR_EXCEPT(true);
    }
  }
}



template <class Scalar> inline 
void Vector<Scalar>::boundscheck(OrdType i, int dim) const
{
  TEST_FOR_EXCEPTION( i < 0 || i >= dim, std::runtime_error,
    "Bounds violation: " << i << "is out of range [0" 
    << ", " << dim << "]");
}


}

#endif
