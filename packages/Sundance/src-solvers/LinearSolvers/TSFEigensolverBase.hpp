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

#ifndef TSFEIGENSOLVERBASE_HPP
#define TSFEIGENSOLVERBASE_HPP

#include "SundanceDefs.hpp"
#include "TSFVectorDecl.hpp" 
#include "TSFSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFLinearOperatorImpl.hpp"

namespace TSFExtended
{
using Teuchos::ParameterList;

/**
 * Base class for eigensolvers for linear eigenvalue problems
 * \f[
 * K x = \lambda M x.
 * \f]
 */
template <class Scalar>
class EigensolverBase 
  : public ObjectWithClassVerbosity<EigensolverBase<Scalar> >
{
public:
  /** */
  EigensolverBase() : params_() {;}

  /** */
  EigensolverBase(const ParameterList& params) : params_(params) {;}

  /** */
  virtual ~EigensolverBase(){;}

  /** 
   * Solve a generalized eigensystem \f$K x = \lambda M x.\f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    const LinearOperator<Scalar>& M,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const = 0 ;

  /** 
   * Solve an eigensystem \f$K x = \lambda x.\f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const 
    {
      LinearOperator<Scalar> M;
      solve(K,M,ev,ew);
    };

  /** 
   * Return the parameter list that was used to define this object. 
   */
  const ParameterList& params() const {return params_;}
  
private:
  ParameterList params_;
};  

}


#endif
