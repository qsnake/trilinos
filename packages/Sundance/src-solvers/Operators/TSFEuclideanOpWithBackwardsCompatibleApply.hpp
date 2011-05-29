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

#ifndef TSFEUCLIDEANOPWITHBACKWARDSCOMPATIBLEAPPLY_HPP
#define TSFEUCLIDEANOPWITHBACKWARDSCOMPATIBLEAPPLY_HPP

#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_DefaultColumnwiseMultiVector.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

/** */
template <class Scalar>
class EuclideanOpWithBackwardsCompatibleApply :
    public virtual SingleScalarEuclideanLinearOpBase<Scalar>
{
public:

  /** */
  EuclideanOpWithBackwardsCompatibleApply(
    const RCP<const ScalarProdVectorSpaceBase<Scalar> >& domain,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> >& range)
    : scalarProdDomain_(domain),
      scalarProdRange_(range)
    {}
    
  /** Implements Euclidean op interface */
  Teuchos::RCP< const ScalarProdVectorSpaceBase<Scalar > > 	
  rangeScalarProdVecSpc () const {return scalarProdRange_;}

  /** Implements Euclidean op interface */
  Teuchos::RCP< const ScalarProdVectorSpaceBase<Scalar > > 	
  domainScalarProdVecSpc () const {return scalarProdDomain_;}

  /** Thyra EuclideanOp apply method */
  void euclideanApply(
    const Thyra::EOpTransp                             transp,
    const Thyra::MultiVectorBase<Scalar>    &X,
    Thyra::MultiVectorBase<Scalar>           *Y,
    const Scalar alpha,
    const Scalar beta
    ) const 
    {
      const Thyra::VectorBase<Scalar>* xVec 
        = dynamic_cast<const Thyra::VectorBase<Scalar>*>(&X);
      
      Thyra::VectorBase<Scalar>* yVec 
        = dynamic_cast<Thyra::VectorBase<Scalar>*>(Y);
      
      if (xVec != 0 && yVec != 0)
      {
        generalApply(transp, *xVec, yVec, alpha, beta);
      }
      else if (xVec == 0 && yVec == 0)
      {
        generalApply(transp, X, Y, alpha, beta);
      }
      else if (X.domain()->dim()==1 && yVec != 0)
      {
        generalApply(transp, *(X.col(0)), yVec, alpha, beta);
      }
      else if (xVec != 0 && Y->domain()->dim()==1)
      {
        generalApply(transp, *xVec, Y->col(0).get(), alpha, beta);
      }
      else
      {
        std::cout << "nX=" << X.domain()->dim()
                  << ", nY=" << Y->domain()->dim() << std::endl;
        TEST_FOR_EXCEPTION(true, std::runtime_error,
          "mix of vectors and multivectors in "
          "OpWithBackwardsCompatibleApply::apply()");
      }
    }

  /** Thyra EuclideanOp apply transpose method */
  void euclideanApplyTranspose(
    const Thyra::EOpTransp                            transp,
    const Thyra::MultiVectorBase<Scalar>    &X,
    Thyra::MultiVectorBase<Scalar>         *Y,
    const Scalar                     alpha,
    const Scalar                     beta
    ) const 
    {
      const Thyra::VectorBase<Scalar>* xVec 
        = dynamic_cast<const Thyra::VectorBase<Scalar>*>(&X);

      Thyra::VectorBase<Scalar>* yVec 
        = dynamic_cast<Thyra::VectorBase<Scalar>*>(Y);

      if (xVec != 0 && yVec != 0)
      {
        generalApply(transp, *xVec, yVec, alpha, beta);
      }
      else if (xVec == 0 && yVec == 0)
      {
        generalApply(transp, X, Y, alpha, beta);
      }
      else if (X.domain()->dim()==1 && yVec != 0)
      {
        generalApply(transp, *(X.col(0)), yVec, alpha, beta);
      }
      else if (xVec != 0 && Y->domain()->dim()==1)
      {
        generalApply(transp, *xVec, Y->col(0).get(), 
          alpha, beta);
      }
      else
      {
        std::cout << "nX=" << X.domain()->dim()
                  << ", nY=" << Y->domain()->dim() << std::endl;
        TEST_FOR_EXCEPTION(true, std::runtime_error,
          "mix of vectors and multivectors in "
          "OpWithBackwardsCompatibleApply::applyTranspose()");
      }
    }

  /** 
   * 
   */
  bool opSupported(Thyra::EOpTransp tr) const 
    {
      return true;
    }


  /**
   * generalApply() applies either the operator or the transpose
   * according to the value of the transpose flag. This method is
   * backwards compatible with TSFCore-based code.
   */
  virtual void generalApply(const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar>    &x,
    Thyra::MultiVectorBase<Scalar>          *y,
    const Scalar            alpha,
    const Scalar            beta) const 
    {
      const Thyra::DefaultColumnwiseMultiVector<Scalar>* xVec 
        = dynamic_cast<const Thyra::DefaultColumnwiseMultiVector<Scalar>*>(&x);
      
      Thyra::DefaultColumnwiseMultiVector<Scalar>* yVec 
        = dynamic_cast<Thyra::DefaultColumnwiseMultiVector<Scalar>*>(y);
      
      TEST_FOR_EXCEPTION(xVec==0, std::runtime_error, 
        "default implementation of "
        "OpWithBackwardsCompatibleApply::generalApply() requires a "
        "DefaultColumnwiseMultiVector");
      
      TEST_FOR_EXCEPTION(yVec==0, std::runtime_error, 
        "default implementation of "
        "OpWithBackwardsCompatibleApply::generalApply() requires a "
        "DefaultColumnwiseMultiVector");
      
      int nXCols = x.domain()->dim();
      int nYCols = y->domain()->dim();
      
      TEST_FOR_EXCEPTION(nXCols != nYCols, std::runtime_error, 
        "mismatched multivector sizes nX=" << nXCols 
        << " and nY=" << nYCols);
      
      for (int i=0; i<nXCols; i++)
      {
        generalApply(M_trans, *(xVec->col(i).get()), 
          (yVec->col(i).get()), alpha, beta);
      }
    }


  

  /**
   *
   */
  virtual void generalApply(const Thyra::EOpTransp M_trans,
    const Thyra::VectorBase<Scalar>    &x,
    Thyra::VectorBase<Scalar>* y,
    const Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
    const Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const = 0 ;
  

private:
  RCP<const ScalarProdVectorSpaceBase<Scalar> > scalarProdDomain_;
  RCP<const ScalarProdVectorSpaceBase<Scalar> > scalarProdRange_;
};

}

#endif
