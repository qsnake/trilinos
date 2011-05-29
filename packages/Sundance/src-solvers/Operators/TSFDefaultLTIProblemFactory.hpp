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

#ifndef TSF_DEFAULT_LTI_PROBLEMFACTORY_HPP
#define TSF_DEFAULT_LTI_PROBLEMFACTORY_HPP


#include "SundanceDefs.hpp"
#include "TSFLTIProblemFactoryBase.hpp"


namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

/** 
 *
 */
template <class Scalar> 
class DefaultLTIProblemFactory
  : public LTIProblemFactoryBase<Scalar>
{
public:
  /** Constructor */
  DefaultLTIProblemFactory(
    int nSteps)
    : LTIProblemFactoryBase<Scalar>(nSteps),
      inputA_(),
      inputC_()
    {}



  /** */
  void init(
    const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& C
    ) 
    {
      inputA_ = A;
      inputC_ = C;
    }

protected:    

  /** Create the operator that advances the system through one timestep. */
  virtual LinearOperator<Scalar> createA() const 
    {
      TEST_FOR_EXCEPT(inputA_.ptr().get() == 0);
      return inputA_;
    }

  /** Create the operator that produces the observable quantities given
   * a state. */
  virtual LinearOperator<Scalar> createC() const 
    {
      TEST_FOR_EXCEPT(inputC_.ptr().get() == 0);
      return inputC_;
    }
  
private:
  LinearOperator<Scalar> inputA_;
  LinearOperator<Scalar> inputC_;
};
}


#endif
