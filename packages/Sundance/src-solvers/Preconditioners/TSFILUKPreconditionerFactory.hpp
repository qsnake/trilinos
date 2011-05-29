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

#ifndef TSFILUKPRECONDITIONERFACTORY_HPP
#define TSFILUKPRECONDITIONERFACTORY_HPP

#include "SundanceDefs.hpp"
#include "TSFPreconditionerFactoryBase.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFILUFactorizableOp.hpp"
#include "TSFLinearSolverBaseDecl.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * 
   */
  template <class Scalar>
  class ILUKPreconditionerFactory
    : public PreconditionerFactoryBase<Scalar>
  {
  public:
    /** Construct with a parameter list */
    ILUKPreconditionerFactory(const ParameterList& params)
      : fillLevels_(1),
        overlapFill_(0),
        relaxationValue_(0.0),
        relativeThreshold_(1.0),
        absoluteThreshold_(0.0),
        leftOrRight_(Right)
    {
      LinearSolverBase<Scalar>::template setParameter<int>(params, &fillLevels_, 
                                                  "Graph Fill");

      LinearSolverBase<Scalar>::template setParameter<int>(params, &overlapFill_, 
                                                  "Overlap");

      LinearSolverBase<Scalar>::template setParameter<double>(params, &relaxationValue_, 
                                                     "Relaxation");

      LinearSolverBase<Scalar>::template setParameter<double>(params, &absoluteThreshold_, 
                                                     "Absolute Threshold");

      LinearSolverBase<Scalar>::template setParameter<double>(params, &relativeThreshold_, 
                                                     "Relative Threshold");

      bool isLeft = false;

      LinearSolverBase<Scalar>::template setParameter<bool>(params, &isLeft, "Left");

      if (isLeft) leftOrRight_ = Left;
      
    }


    /** virtual dtor */
    virtual ~ILUKPreconditionerFactory(){;}

    
    /** */
    virtual Preconditioner <Scalar>
    createPreconditioner(const LinearOperator<Scalar>& A) const 
    {
      /* In order for ILU factorization to work, the operator A must
       * implement the ILUFactorizableOp interface. We cast A's pointer
       * to a ILUFactorizableOp ptr. If the cast fails, throw a spoke. */
      
      const ILUFactorizableOp<Scalar>* fop 
        = dynamic_cast<const ILUFactorizableOp<Scalar>*>(A.ptr().get());

      TEST_FOR_EXCEPTION(fop==0, std::runtime_error,
                         "ILUKPreconditionerFactory attempted to "
                         "create an ILU preconditioner for an operator type "
                         "that does not implement the ILUFactorizableOp "
                         "interface. The op is " << A.description());

      
      /* Now we can delegate the construction of the ILU factors to 
      * the factorizable op. */
      Preconditioner<Scalar> P;
      fop->getILUKPreconditioner(fillLevels_,
                                 overlapFill_,
                                 relaxationValue_,
                                 relativeThreshold_,
                                 absoluteThreshold_,
                                 leftOrRight_,
                                 P);
      /* Return the preconditioner */
      return P;
    }

    /* Handleable boilerplate */
    GET_RCP(PreconditionerFactoryBase<Scalar>);
  private:

    int fillLevels_;
    int overlapFill_;
    Scalar relaxationValue_;
    Scalar relativeThreshold_;
    Scalar absoluteThreshold_;
    LeftOrRight leftOrRight_;
  };


}

#endif
