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

#ifndef TSFNONLINEAROPERATOR_HPP
#define TSFNONLINEAROPERATOR_HPP

#include "SundanceDefs.hpp"
#include "SundanceHandle.hpp"
#include "TSFNonlinearOperatorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace TSFExtended
{
  using Thyra::Index;
  using namespace Teuchos;

  /** 
   * User-level nonlinear operator class
   */
  template <class Scalar>
  class NonlinearOperator : public Sundance::Handle<NonlinearOperatorBase<Scalar> >
    {
    public:
      /* boilerplate ctors */
      HANDLE_CTORS(NonlinearOperator<Scalar>, NonlinearOperatorBase<Scalar>);

      /** */
      VectorSpace<Scalar> domain() const 
      {return this->ptr()->domain();}

      /** */
      VectorSpace<Scalar>  range() const 
      {return this->ptr()->range();}

      /** */
      void setEvalPt(const Vector<double>& evalPt)
      {
        this->ptr()->setEvalPt(evalPt);
      }
      
      /** */
      LinearOperator<Scalar> getJacobian() const 
      {
        return this->ptr()->getJacobian();
      }

      /** */
      Vector<double> getFunctionValue() const 
      {
        return this->ptr()->getFunctionValue();
      }

      

      /** */
      Vector<double> getInitialGuess() const 
      {
        return this->ptr()->getInitialGuess();
      }

      /** */
      Vector<double> currentEvalPt() const 
      {
        return this->ptr()->currentEvalPt();
      }

    private:
    };
}


#endif
