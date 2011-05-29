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

#ifndef TSFBCPARTITIONEDVSBUILDER_HPP
#define TSFBCPARTITIONEDVSBUILDER_HPP

#include "TSFVectorSpaceImpl.hpp"
#include "TSFVectorType.hpp"
#include "Teuchos_Array.hpp"
#include "TSFProductVectorSpaceImpl.hpp"

namespace TSFExtended
{
using namespace Teuchos;

/**
 * This helper function builds a product space in which BC and internal
 * rows are partitioned into physically distinct spaces
 */
template <class Scalar>
  VectorSpace<Scalar> buildPartitionedSpace(
    int nTotalDofs,
    int lowestLocalDof,
    int nLocalDofs,
    const Array<int>& isBCIndex,
    const VectorType<Scalar>& internalType,
    const VectorType<Scalar>& bcType,
    const MPIComm& comm
    )
  {
    int nBCDofs = 0;
    for (int i=0; i<nLocalDofs; i++)
    {
      if (isBCIndex[i]) nBCDofs++;
    }

    /* sum number of BC Dofs over all processors */
    int nTotalBCDofs = nBCDofs;
    comm.allReduce(&nBCDofs, &nTotalBCDofs, 1, MPIComm::INT, MPIComm::SUM);
    int nTotalInteriorDofs = nTotalDofs - nTotalBCDofs;


    Array<int> interiorDofs(nLocalDofs - nBCDofs);
    Array<int> bcDofs(nBCDofs);
    int iBC = 0;
    int iIn = 0;

    for (int i=0; i<nLocalDofs; i++)
    {
      if (isBCIndex[i]) bcDofs[iBC++] = lowestLocalDof+i;
      else interiorDofs[iIn++] = lowestLocalDof+i;
    }

    int p = MPIComm::world().getRank();


    VectorSpace<double> bcSpace = bcType.createSpace(nTotalBCDofs, nBCDofs,
      &(bcDofs[0]), comm);
    VectorSpace<double> interiorSpace = internalType.createSpace(nTotalInteriorDofs, nLocalDofs-nBCDofs,
      &(interiorDofs[0]), comm);

    return productSpace<double>(interiorSpace, bcSpace);
  }

}


#endif
