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

#ifndef TSFANASAZIEIGENSOLVER_DECL_HPP
#define TSFANASAZIEIGENSOLVER_DECL_HPP

#include "SundanceDefs.hpp"
#include "TSFVectorDecl.hpp" 
#include "TSFSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFEigensolverBase.hpp"

namespace TSFExtended
{
using Teuchos::ParameterList;

/**
 * Object wrapper for Anasazi eigenvalue solver.
 */
template <class Scalar>
class AnasaziEigensolver
  : public EigensolverBase<Scalar>,
    public Sundance::Handleable<EigensolverBase<Scalar> >
{
public:
  /** */
  AnasaziEigensolver(const ParameterList& params) 
    : EigensolverBase<Scalar>(params) {;}

  /**
   * Solve a generalized eigenvalue problem \f$ K x = \lambda M x \f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    const LinearOperator<Scalar>& M,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const ;

  /** \name Handleable interface */
  //@{
  /** Return a ref counted pointer to a newly created object */
  virtual RCP<EigensolverBase<Scalar> > getRcp() 
    {return rcp(this);}
  //@}
  
private:

  static Time& solveTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("AnasaziEigensolver solve()"); 
      return *rtn;
    }

  static Time& precondBuildTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("AnasaziEigensolver building preconditioner"); 
      return *rtn;
    }

};




}


#endif
