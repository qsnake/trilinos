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

#ifndef TSFVECTORSPACEIMPL_HPP
#define TSFVECTORSPACEIMPL_HPP


#include "Thyra_ProductVectorSpaceBase.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFSequentialIteratorImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;


static inline Time& createVecTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("vector allocation"); 
  return *rtn;
}

 
//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::operator==(const VectorSpace<Scalar>& other) const 
{
  return isCompatible(other);  
}


//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::operator!=(const VectorSpace<Scalar>& other) const 
{
  return !(operator==(other));
}
    


//========================================================================
template <class Scalar>
Vector<Scalar> VectorSpace<Scalar>::createMember() const 
{
  TimeMonitor timer(createVecTimer());
  Vector<Scalar> rtn = Thyra::createMember(this->ptr());
  rtn.setToConstant(0.0);
  return rtn;
}
    


//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::lowestLocallyOwnedIndex() const
{
  const Thyra::SpmdVectorSpaceBase<Scalar>* mpiSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
  if (mpiSpace != 0)
    {
      return mpiSpace->localOffset();
    }
  const Thyra::SpmdVectorSpaceBase<Scalar>* serialSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
   if (serialSpace != 0)
     {
       return 0;
     }
   TEST_FOR_EXCEPTION(mpiSpace == 0 && serialSpace==0, std::runtime_error,
		      "don't know how to compute lowest local index for "
		      "a vector space that is neither MPI nor serial");
   return 0;
}

//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::numLocalElements() const
{
  if (numBlocks() > 1)
  {
    int rtn = 0;
    for (int b=0; b<numBlocks(); b++) 
    {
      rtn += getBlock(b).numLocalElements();
    } 
    return rtn;
  }
  
  const Thyra::SpmdVectorSpaceBase<Scalar>* mpiSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
  if (mpiSpace != 0)
    {
      return mpiSpace->localSubDim();
    }
   const Thyra::SpmdVectorSpaceBase<Scalar>* serialSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
   if (serialSpace != 0)
     {
       return dim();
     }
   TEST_FOR_EXCEPTION(mpiSpace == 0 && serialSpace==0, std::runtime_error,
		      "don't know how to compute number of local elements for "
		      "a vector space that is neither MPI nor serial");
   return 0;
}
    



//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::isCompatible(const VectorSpace<Scalar>& vecSpc) const 
{
  TEST_FOR_EXCEPTION(vecSpc.ptr().get() == 0, std::runtime_error,
                     "null argument in VectorSpace<Scalar>::isCompatible()");
  return this->ptr().get()->isCompatible(*(vecSpc.ptr().get()));
}





//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::contains(const Vector<Scalar> &vec) const
{
  return (operator==(vec.space()));
}


//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::numBlocks() const
{
  const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
    dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->ptr().get());
  if (pvs != 0)
    {
      return pvs->numBlocks();
    }
  return 1;
}



//========================================================================
template <class Scalar>
VectorSpace<Scalar> VectorSpace<Scalar>::getBlock(const int i) const
{
  const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
    dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->ptr().get());
  TEST_FOR_EXCEPTION(pvs == 0 && numBlocks()!=1, std::runtime_error,
		     "Space not a ProductVectorSpace" << std::endl);
  if (pvs != 0)
    {
      return pvs->getBlock(i);
    }
  return *this;
}


// //========================================================================
// template <class Scalar>
// void VectorSpace<Scalar>::setBlock(int i, 
// 				   const VectorSpace<Scalar>& space)
// {
//   const Thyra::ProductVectorSpace<Scalar>*  pvs = 
//     dynamic_cast<const Thyra::ProductVectorSpace<Scalar>* >  (this->ptr().get());

//   TEST_FOR_EXCEPTION(pvs == 0, std::runtime_error,
// 		     "Can't set block of vector space that is " <<
// 		     "not a ProductVectorSpace.");

//   Thyra::ProductVectorSpace<Scalar>* pvsc = const_cast<ProductVectorSpace<Scalar>*> (pvs);
//   pvsc->setBlock(i, space);
// }


//========================================================================
template <class Scalar> inline
SequentialIterator<Scalar> VectorSpace<Scalar>::begin() const
{
  OrdType blockIndex = -1;
  OrdType indexInCurrentBlock = -1;
  OrdType globalIndex = -1;

  /* we need to check for remaining data to deal with the case where 
   * a space is empty */
  bool dataRemains = advanceIndex(blockIndex, indexInCurrentBlock, globalIndex);

  if (dataRemains)
  {
    return SequentialIterator<Scalar>(*this, blockIndex, indexInCurrentBlock, globalIndex);
  }
  else
  {
    return end();
  }
}


//========================================================================
template <class Scalar> inline
SequentialIterator<Scalar> VectorSpace<Scalar>::end() const
{
  return SequentialIterator<Scalar>(*this);
}





template <class Scalar> inline
bool VectorSpace<Scalar>::advanceIndex(
  OrdType& blockIndex, 
  OrdType& indexInCurrentBlock,
  OrdType& globalIndex) const 
{
  /* block index == -1 indicates the initialization call */
  if (blockIndex < 0)
  {
    /* Find the start of the first nonempty block */
    for (int i=0; i<numBlocks(); i++)
    {
      if (getBlock(i).dim()==0) continue;
      blockIndex = i;
      indexInCurrentBlock = 0;
      globalIndex = lowestLocallyOwnedIndex();
      return true;
    }
    /* if we've made it to this point, all blocks are empty so there's
     * no data. */
    return false;
  }

  /* If we're not a product space, then all we need to do is increment
   * the index and then check whether we've run off the end. 
   * Note: at this point, we can only deal with unit stride increments.
   * To do anything else would require more intimate communication with
   * the underlying VectorSpaceBase object, perhaps through something
   * like an advanceIndex() function. 
   */
  bool isProductSpace = (0 != 
    dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->ptr().get()));
  if (!isProductSpace)
  {
    indexInCurrentBlock++;
    globalIndex++;
    if (indexInCurrentBlock >= numLocalElements())
    {
      return false;
    }
    return true;
  }
  else 
  {
    /* If we are a product space, first try to advance within 
     * the current block. */
    OrdType subBlock = 0;
    if (getBlock(blockIndex).advanceIndex(subBlock, indexInCurrentBlock, globalIndex))
    {
      /* Advance was successful. */
      return true;
    }
    else
    {
      /* If no data remains in the current block, find the next block
       * the contains data, and restart the indexing within that block. 
       */
      for (int i=blockIndex+1; i<numBlocks(); i++)
      {
        if (getBlock(i).dim()==0) continue;
        indexInCurrentBlock = 0;
        globalIndex++;
        blockIndex = i;
        return true;
      }
      /* no data remaining in any block */
      return false;
    }
  }
}







#endif
