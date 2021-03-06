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

#ifndef LOCA_MULTICONTINUATION_EXTENDEDVECTOR_H
#define LOCA_MULTICONTINUATION_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    /*! 
      \brief %Vector class to hold solution vectors, %Newton vectors, etc.
      for continuation equations.
    */
    /*! 
       This class uses the LOCA::Extended::Vector implementation to store the
       solution and parameter components of the continuation vector
       and merely provides an interface for naming which components of the
       multivector these quantities correspond to.
    */
    class ExtendedVector : public LOCA::Extended::Vector {

      /*! 
       * \brief Declare LOCA::MultiContinuation::ExtendedMultiVector as a 
       * friend class so it can call protected methods.
       */
      friend class ExtendedMultiVector;

    public:
  
      //! Constructor.
      ExtendedVector(const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     const NOX::Abstract::Vector& xVec, int nScalars);

      //! Copy constructor.
      ExtendedVector(const ExtendedVector& source, 
		     NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ExtendedVector();

      //! Assignment operator
      virtual ExtendedVector& 
      operator=(const ExtendedVector& y);

      //! Assignment operator
      virtual LOCA::Extended::Vector& 
      operator=(const LOCA::Extended::Vector& y);

      //! Assignment operator
      virtual NOX::Abstract::Vector& 
      operator=(const NOX::Abstract::Vector& y);

      //! Assignment operator
      virtual Teuchos::RCP<NOX::Abstract::Vector>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Returns the solution vector component of extended vector
      virtual Teuchos::RCP<const NOX::Abstract::Vector> 
      getXVec() const;

      //! Returns the solution vector component of extended vector
      virtual Teuchos::RCP<NOX::Abstract::Vector> getXVec();

    protected:

      //! Constructor
      ExtendedVector(const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     int nScalars);

      //! generate an extended multi-vector
      /*!
       * Returns an empty multivector of type 
       * LOCA::MultiContinuation::ExtendedMultiVector
       */
      virtual Teuchos::RCP<LOCA::Extended::MultiVector>
      generateMultiVector(int nColumns, int nVectorRows, 
			  int nScalarRows) const;

    }; // class ExtendedVector
  } // namespace Continuation
} // namespace LOCA

#endif
