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

#ifndef TSF_SERIAL_GHOSTVIEW_HPP
#define TSF_SERIAL_GHOSTVIEW_HPP

#include "SundanceDefs.hpp"
#include "TSFGhostView.hpp"
#include "TSFSerialVector.hpp"
#include "Teuchos_Utils.hpp"



namespace TSFExtended
{
  using namespace Teuchos;


  /**
   * Dummy ghost element viewer for serial vectors. 
   */
  class SerialGhostView : public GhostView<double>
    {
    public:
      /** */
      SerialGhostView(const RCP<SerialVector>& vec)
        : vec_(vec) 
      {;}

      /** virtual dtor */
      virtual ~SerialGhostView(){;}

      /** Indicate whether the given global index is accessible in this view */
      bool isAccessible(OrdType globalIndex) const 
      {return true;}

      /** get the element at the given global index */
      const double& getElement(OrdType globalIndex) const 
        {return vec_->getElement(globalIndex);}

      /** get the batch of elements at the given global indices */
      void getElements(const OrdType* globalIndices, OrdType numElems,
                       Array<double>& elems) const 
        {
          vec_->getElements(globalIndices, numElems, elems);
        }

      /** */
      void print(std::ostream& os) const {vec_->print(os);}
    private:
      RCP<const SerialVector> vec_;
    };
  
}

#endif
