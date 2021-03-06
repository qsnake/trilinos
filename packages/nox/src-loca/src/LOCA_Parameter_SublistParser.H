// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_PARAMETER_SUBLISTPARSER_H
#define LOCA_PARAMETER_SUBLISTPARSER_H

#include "Teuchos_RCP.hpp"

#include "NOX_Common.H"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Parameter {

    //! Class to parse a parameter list for sublists
    /*!
     * This class parses a supplied parameter list and locates various
     * sublists.  This saves the code from having to traverse the parameter
     * list to find sublists itself, and puts in one location the
     * hard-coded structure of the parameter list.
     */
    class SublistParser {

      //! Map type to hold sublists
      typedef map<string, 
		  Teuchos::RCP<Teuchos::ParameterList> > SublistMap;

      //! Iterator type to traverse SublistMap
      typedef SublistMap::iterator SublistMapIterator;

    public:

      //! Constructor
      SublistParser(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      ~SublistParser();

      //! Parse parameter list to find sublists
      void parseSublists(
	     const Teuchos::RCP<Teuchos::ParameterList>& topLevelParams);

      //! Return sublist of name \c name
      Teuchos::RCP<Teuchos::ParameterList> 
      getSublist(const string& name);

    private:

      //! Private to prohibit copying
      SublistParser(const SublistParser&);

      //! Private to prohibit copying
      SublistParser& operator = (const SublistParser&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Map to hold sublists
      SublistMap sublistMap;

    }; // Class SublistParser

  } // Namespace Parameter

} // Namespace LOCA

#endif
