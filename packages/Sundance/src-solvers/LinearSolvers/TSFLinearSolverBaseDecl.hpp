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

#ifndef TSFLINEARSOLVERBASEDECL_HPP
#define TSFLINEARSOLVERBASEDECL_HPP

#include "SundanceDefs.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceNamedObject.hpp"
#include "TSFSolverState.hpp"
#include "Teuchos_ParameterList.hpp"

namespace TSFExtended
{
  using namespace Teuchos;
  template <class Scalar>
  class LinearOperator;

  template <class Scalar>
  class Preconditioner;

  template <class Scalar>
  class PreconditionerFactory;

  template <class Scalar>
  class Vector;
  

  /** */
  template <class Scalar>
  class LinearSolverBase : public DefaultObjectWithVerbosity,
                           public NamedObject
  {
  public:
    /** */
    LinearSolverBase(const ParameterList& params);

    /** */
    virtual ~LinearSolverBase(){;}

    /** */
    virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                                      const Vector<Scalar>& rhs,
                                      Vector<Scalar>& soln) const = 0;

    /** Change the convergence tolerance. Default does nothing. */
    virtual void updateTolerance(const double& tol) {;}

    /** Set a user-defined preconditioning operator. Default is an error. */
    virtual void setUserPrec(const PreconditionerFactory<Scalar>& pf);

    /** Set a user-defined preconditioning operator. Default is an error. */
    virtual void setUserPrec(const LinearOperator<Scalar>& P,
      const LinearSolver<Scalar>& pSolver);

    /** */
    const ParameterList& parameters() const ;

    /** */
    ParameterList& parameters();

    /** */
    static std::string verbosityParam();

    /** */
    template <typename T>
    static void setParameter(const ParameterList& params,
                             T* valuePtr, 
                             const std::string& paramName);
  private:
    ParameterList params_;
  };
}

#endif
