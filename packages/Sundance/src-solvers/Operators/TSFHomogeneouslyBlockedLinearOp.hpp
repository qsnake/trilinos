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

#ifndef TSFHOMOGENEOUSLYBLOCKEDLINEAROP_HPP
#define TSFHOMOGENEOUSLYBLOCKEDLINEAROP_HPP

#include "SundanceDefs.hpp"
#include "TSFSimplifiedLinearOpBase.hpp"

namespace TSFExtended
{
using namespace Teuchos;

/**
 * HomogeneouslyBlockedLinearOp is a helper class providing a convenient
 * way to build operators having a block structure in which all blocks
 * have the same domain and range. Such structures arise in, for example,
 * space-time operators or stochastic projection methods.
 *
 * @author Kevin Long (krlong@sandia.gov)
 */
template <class Scalar>
class HomogeneouslyBlockedLinearOp :
    public virtual LinearOpBase<Scalar>
{
public:

  /*
   * Construct by specifying the domain and range spaces for a single
   * block and the number of blocks over which to replicate this.
   */
  HomogeneouslyBlockedLinearOp(
    const VectorSpace<Scalar>& singleBlockDomain, 
    int numDomainBlocks,
    const VectorSpace<Scalar>& singleBlockRange,
    int numRangeBlocks
    )
    : singleBlockDomain_(singleBlockDomain),
      numDomainBlocks_(numDomainBlocks),
      singleBlockRange_(singleBlockRange),
      numRangeBlocks_(numRangeBlocks),
      domain_(),
      range_()
    {
      Array<VectorSpace<Scalar> > d(numDomainBlocks_, singleBlockDomain_);
      Array<VectorSpace<Scalar> > r(numRangeBlocks_, singleBlockRange_);
      domain_ = productSpace(d);
      range_ = productSpace(r);
    }      

  /** 
   * \brief Return a smart pointer for the range space 
   * for <tt>this</tt> operator.
   */
  Teuchos::RCP< const Thyra::VectorSpaceBase<Scalar> > range() const 
    {
      return range_;
    }

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   */
  Teuchos::RCP< const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {
      return domain_;
    }

protected:

  /** Return the number of block rows */
  int numBlockRows() const {return numRangeBlocks_;}

  /** Return the number of block columns */
  int numBlockCols() const {return numDomainBlocks_;}

  /** Return the domain space of a single block. By construction
   * this is the same for all blocks.  */
  const VectorSpace<Scalar>& singleBlockDomain() const 
    {
      return singleBlockDomain_;
    }

  /** Return the range space of a single block. By construction
   * this is the same for all blocks.  */
  const VectorSpace<Scalar>& singleBlockRange() const 
    {
      return singleBlockRange_;
    }

private:
  VectorSpace<Scalar> singleBlockDomain_;
  int numDomainBlocks_;
  VectorSpace<Scalar> singleBlockRange_;
  int numRangeBlocks_;
  Teuchos::RCP< const Thyra::VectorSpaceBase<Scalar> > domain_;
  Teuchos::RCP< const Thyra::VectorSpaceBase<Scalar> > range_;
}; 
}

#endif
