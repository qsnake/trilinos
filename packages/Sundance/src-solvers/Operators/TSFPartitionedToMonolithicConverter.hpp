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

#ifndef TSFPARTITIONED_TO_MONOLITHIC_CONVERTER_HPP
#define TSFPARTITIONED_TO_MONOLITHIC_CONVERTER_HPP

#include "TSFVectorSpaceImpl.hpp"
#include "TSFVectorImpl.hpp"
#include "TSFVectorType.hpp"
#include <set>

namespace TSFExtended
{
  using namespace Teuchos;
  using namespace Thyra;

  /** */
  class PartitionedToMonolithicConverter
  {
  public:

    /**  */
    PartitionedToMonolithicConverter(
      const VectorSpace<double>& partitionedSpace,
      const RCP<Array<int> >& isBCCol,
      const VectorSpace<double>& monolithicSpace
      ) : partitionedSpace_(partitionedSpace),
          monolithicSpace_(monolithicSpace),
          indices_(2)
      {
        int n = monolithicSpace.numLocalElements();
        int low = monolithicSpace.lowestLocallyOwnedIndex();

        indices_[0].reserve(n);
        indices_[1].reserve(n);

        const Array<int>& bc = *isBCCol;

        for (int i=0; i<n; i++)
        {
          indices_[bc[i]].append(low+i);
        }
      }

    /** */
    void convert(
      const Vector<double>& in,
      Vector<double>& out) const 
      {
        TEST_FOR_EXCEPTION(in.space().numBlocks() != 2, std::runtime_error,
          "PartitionedToMonolithicConverter::convert() expects 2 blocks "
          "in input vector. Input vector has numBlocks()="
          << in.space().numBlocks());
        Array<Vector<double> > x(2);
        x[0] = in.getBlock(0);
        x[1] = in.getBlock(1);

        cout << "index array=" << indices_ << std::endl;

        for (int b=0; b<2; b++)
        {
          VectorSpace<double> space = x[b].space();
          int j=0;
          for (SequentialIterator<double> i=space.begin(); i!=space.end(); i++,j++)
          {
            double val = x[b][i];
            int globalIndex = indices_[b][j];
            out.setElement(globalIndex, val);
          }
        }
      }

  protected:

  private:
    VectorSpace<double> partitionedSpace_;
    VectorSpace<double> monolithicSpace_;
    Array<Array<int> > indices_;
  };
}


#endif
