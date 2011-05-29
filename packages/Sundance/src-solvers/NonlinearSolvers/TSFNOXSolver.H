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

#ifndef TSFNOXSOLVER_HPP
#define TSFNOXSOLVER_HPP

#include "SundanceDefs.hpp"
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"
#include "NOX_TSF_StatusTestBuilder.H"
#include "NOX_Multiphysics_Solver_Manager.H"
#include "Teuchos_TestForException.hpp"   
#include "Teuchos_ParameterList.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   *
   */
  class NOXSolver 
  {
  public:
    /** */
    NOXSolver(){;}
    /** */
    NOXSolver(const ParameterList& params);
    /** */
    NOXSolver(const ParameterList& nonlinParams,
      const LinearSolver<double>& linSolver);

    /** */
    NOX::StatusTest::StatusType solve(const NonlinearOperator<double>& F, 
      Vector<double>& soln) const ;

    /** */
    const LinearSolver<double>& linSolver() const 
      {return linSolver_;}


  private:

    LinearSolver<double> linSolver_;
    mutable RCP<NOX::StatusTest::Generic> statusTest_;
    mutable ParameterList params_;
    mutable ParameterList printParams_;
  };
}

#endif
