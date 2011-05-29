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

#ifndef TSFEPETRAGHOSTVIEW_HPP
#define TSFEPETRAGHOSTVIEW_HPP

#include "SundanceDefs.hpp"
#include "TSFGhostImporter.hpp"
#include "TSFGhostView.hpp"
#include "Epetra_Vector.h"
#include "Teuchos_Utils.hpp"



namespace TSFExtended
{
  using namespace Teuchos;


  /**
   * Ghost element viewer for Epetra vectors
   */
  class EpetraGhostView : public GhostView<double>
    {
    public:
      /** */
      EpetraGhostView()
        : ghostView_() 
      {;}

      /** virtual dtor */
      virtual ~EpetraGhostView(){;}

      /** Indicate whether the given global index is accessible in this view */
      bool isAccessible(OrdType globalIndex) const 
      {return ghostView_->Map().MyGID(globalIndex);}

      /** get the element at the given global index */
      const double& getElement(OrdType globalIndex) const ;

      /** get the batch of elements at the given global indices */
      void getElements(const OrdType* globalIndices, int numElems,
                       Array<double>& elems) const ;

      /** */
      void import(const Epetra_Import& importer,
                  const Epetra_Vector& srcObject);

      /** */
      void print(std::ostream& os) const ;
    private:
      RCP<Epetra_Vector> ghostView_;
    };
  
}

#endif
