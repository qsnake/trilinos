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

#ifndef TSFPRECONDITIONER_HPP
#define TSFPRECONDITIONER_HPP

#include "SundanceDefs.hpp"
#include "SundanceHandle.hpp"
#include "TSFPreconditionerBase.hpp"

namespace TSFExtended
{
  /**
   * 
   */
  template <class Scalar> 
  class Preconditioner : public Sundance::Handle<PreconditionerBase<Scalar> >
  {
  public:
    /* Boilerplate ctors */
    HANDLE_CTORS(Preconditioner, PreconditionerBase<Scalar>);

    /** Change the value of a double parameter */
    void changeParameter(const std::string& name, const double& value);

    /** Change the value of an integer parameter */
    void changeParameter(const std::string& name, int value);

    
    
    /** Left preconditioner */
    LinearOperator<Scalar> left() const ;
    
    /** Right preconditioner */
    LinearOperator<Scalar> right() const ;
    
    /** return true if this preconditioner has both left and
     * right components. */
    bool isTwoSided() const {return hasLeft() && hasRight();}
    
    /** return true if this preconditioner has a nontrivial left component */
    bool hasLeft() const ;
    
    /** return true if this preconditioner has
     * a nontrivial right component */
    bool hasRight() const ;
    
    /** return true if this preconditioner has neither left nor
     * right operators defined */
    bool isIdentity() const {return !hasLeft() && !hasRight();}
  };

  

  template <class Scalar> inline 
  LinearOperator<Scalar> Preconditioner<Scalar>::left() const 
  {
    TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Preconditioner<Scalar>::left()");
    return this->ptr()->left();
  }

  template <class Scalar> inline 
  LinearOperator<Scalar> Preconditioner<Scalar>::right() const 
  {
    TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Preconditioner<Scalar>::right()");
    return this->ptr()->right();
  }

  template <class Scalar> inline
  bool Preconditioner<Scalar>::hasLeft() const 
  {
    return (this->ptr().get()!=0 && this->ptr()->hasLeft());
  }

  template <class Scalar> inline
  bool Preconditioner<Scalar>::hasRight() const 
  {
    return (this->ptr().get()!=0 && this->ptr()->hasRight());
  }

  
}

#endif
