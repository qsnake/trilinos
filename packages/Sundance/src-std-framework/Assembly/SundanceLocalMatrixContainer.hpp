/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_LOCALMATRIXCONTAINER_H
#define SUNDANCE_LOCALMATRIXCONTAINER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBuilder.hpp"

namespace Sundance
{
/** 
 * LocalMatrixContainer is a container for local matrix and vector values.
 * Local matrices and vectors are stored in identical ways: as 1D arrays
 * of doubles. toi
 *
 * Each integrator gets its own LocalMatrixContainer object, with information on
 * term grouping specific to the terms handled by that integrator. To conserve
 * memory, all LocalMatrixContainer objects share a common pool of data vectors. After 
 * insertion, the common pool can be reused in the next integration. It is the 
 * responsibility of the integrator to zero out and resize each data vector
 * as needed.
 *   
 * To leave the system as flexible as possible, this class does not specify
 * anything about the ordering of elements in the matrix. It is the
 * responsibility of the integrator and inserter classes to use a consistent
 * ordering. Different integrator-inserter combinations may well use different
 * orderings.
 *
 * It is often the case that several terms in a PDE will yield local 
 * matrices that have the same values, perhaps differing by a constant, 
 * but which are associated with different test and unknown functions. For this reason,
 * each batch of local matrix values has associated with it arrays
 * of (test, unk) pairs and coefficients.  
 */
class LocalMatrixContainer
{
public:
  /** */
  LocalMatrixContainer(const Array<int>& isTwoForm,
    const Array<Array<int> >& testID,
    const Array<Array<int> >& unkID,
    const Array<Array<double> >& coeffs);

  /** Return the data vector for the i-th batch */
  const RCP<Array<double> >& dataVector(int i) const 
    {return workspace()[i];}

  /** Indicate whether the i-th batch is a two form */
  bool isTwoForm(int i) const {return isTwoForm_[i];}

  /** Return the array of testIDs whose local matrix values are grouped int
   * the i-th batch */
  const Array<int>& testID(int i) const {return testID_[i];}

  /** Return the array of unkIDs whose local matrix values are grouped int
   * the i-th batch */
  const Array<int>& unkID(int i) const {return unkID_[i];}

  /** Return the array of coefficients to be used with the i-th batch */
  const Array<double>& coeffs(int i) const {return coeffs_[i];}

private:
      
  /** */
  static Array<RCP<Array<double> > >& workspace() 
    {static Array<RCP<Array<double> > > rtn; return rtn;}

  /** */
  Array<int> isTwoForm_;

  /** */
  Array<Array<int> > testID_;

  /** */
  Array<Array<int> > unkID_;

  /** */
  Array<Array<double> > coeffs_;
};
}



#endif
