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

#ifndef TSF_SERIAL_GHOSTIMPORTER_HPP
#define TSF_SERIAL_GHOSTIMPORTER_HPP

#include "SundanceDefs.hpp"
#include "TSFGhostImporter.hpp"
#include "TSFSerialGhostView.hpp"
#include "Teuchos_Utils.hpp"



namespace TSFExtended
{
  using namespace Teuchos;


  /**
   * Ghost element importer for serial vectors. This class doesn't have
   * much to do, but is necessary to maintain a consistent interface.
   */
  class SerialGhostImporter : public GhostImporter<double>
    {
    public:
      /** */
      SerialGhostImporter(){;}
      /** virtual dtor */
      virtual ~SerialGhostImporter() {;}

      /** 
       * Import the ghost elements of the given vector
       * as specified during construction of this object. 
       */
      virtual void importView(const Vector<double>& x,
                              RCP<GhostView<double> >& ghostView) const ;

    private:
      
    };
  
}

#endif
