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
                                                                                
//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NLS_PetraGroup to provide for 
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>
#include "petscksp.h"
#include "NOX_Petsc_Interface.H"

// ---------- Forward Declarations ----------
class FiniteDifference;

class  Problem_Interface : public NOX::Petsc::Interface
{
public:
  Problem_Interface(FiniteDifference& Problem);
  ~Problem_Interface();

  //! Compute and return RHS.  Returns true if computation was successful.
  bool computeF(const Vec& x, Vec& RHS);

  //! Compute an explicit Jacobian.  Returns true if computation was successful.
  bool computeJacobian(const Vec& x, Mat& Jac);

  //! Compute the matrix M that will be used as the preconditioner.  Returns true if computation was successful.
  bool computePreconditioner(Mat& M);

  //! Return the action of the preconditioner on a vector.  Returns true if computation was successful.
  bool preconditionVector(Vec& y);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 
  FiniteDifference& problem;
};

#endif

