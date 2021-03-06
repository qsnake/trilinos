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

#ifndef LOCA_LAPACK_FACTORY_H
#define LOCA_LAPACK_FACTORY_H

#include "LOCA_Abstract_Factory.H"    // base class

// Forward declarations
namespace LOCA{
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace LAPACK {

    //! Implementation of the LOCA::Abstract::Factory for LAPACK groups.
    class Factory : public LOCA::Abstract::Factory {

    public:

      //! Constructor
      Factory();

      //! Destructor
      virtual ~Factory();

      //! Initialize factory
      virtual void 
      init(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      /*! 
       * @name Strategy create methods
       */
      //@{

       //! Create bordered system solver strategy
      virtual bool
      createBorderedSolverStrategy(
       const string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy);

      //! Create eigensolver strategy
      virtual bool
      createEigensolverStrategy(
         const string& strategyName,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	 Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>& strategy);

      //@}

    private:

      //! Private to prohibit copying
      Factory(const Factory& fac);

      //! Private to prohibit copying
      Factory& operator = (const Factory& fac);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace LAPACK

} // Namespace LOCA

#endif
