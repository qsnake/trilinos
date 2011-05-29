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

#ifndef TSFACCESSIBLEVECTOR_HPP
#define TSFACCESSIBLEVECTOR_HPP

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"


namespace TSFExtended
{
  

  /**
   * TSFExtended::AccessibleVector defines an interface through which
   * elements for a vector can be accessed. Element access is occasionally
   * used by application codes in probing results vectors, 
   * but should rarely be used by high-performance solver codes; this 
   * capability is therefore in TSFExtended rather than Thyra.
   *
   * @author Kevin Long (krlong@sandia.gov)
   */
  template <class Scalar>
  class AccessibleVector 
    {
    public:
      /** virtual dtor */
      virtual ~AccessibleVector() {;}

      /** get the element at the given global index */
      virtual const Scalar& getElement(OrdType globalIndex) const = 0 ;

      /** get a batch of elements. Slow default implementation loops
       * over calls to getElement(). */

      virtual void getElements(const int* globalIndices, int numElems,
        Teuchos::Array<Scalar>& elems) const 
        {
          elems.resize(numElems);
          for (int i=0; i<numElems; i++)
            {
              elems[i] = getElement(globalIndices[i]);
            }
        }
    };
}



#endif
