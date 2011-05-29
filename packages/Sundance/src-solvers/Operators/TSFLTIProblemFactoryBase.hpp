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

#ifndef TSFLTIPROBLEMFACTORYBASE_HPP
#define TSFLTIPROBLEMFACTORYBASE_HPP

#include "SundanceDefs.hpp"
#include "TSFInverseLTIOp.hpp"


namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

/** 
 * LTIProblemFactoryBase is an interface for building 
 * operators related to reduced-order
 * LTI problems. The operators are
 *
 * <UL>
 * <LI> \f$ A \f$ -- operator defining a single step of 
 * the discrete time dynamical system
 * <LI> \f$ C \f$ -- operator defining the transformation from state
 * variables to observed variables at a single timestep
 * <LI> \f$ {\hat A}^{-1} \f$ -- block operator that produces a sequence
 * of state variables. This will be of type InverseLTIOp
 * <LI> \f$ {\hat C} \f$ -- block diagonal consisting of observation operators
 * for all timesteps.
 * <LI> \f$ F \f$ -- operator that maps the initial state to a full
 * space-time state vector padded with zeros: 
 * \f$ F = \left[\begin{array}{cccc} I & 0 & 0 & \cdots\end{array}\right]\f$
 * </UL>
 *
 * <H1> Create() methods versus get() methods </H1>
 * 
 * The create methods actually construct operators. They are protected and
 * can therefore only be called by internal methods. Client access to
 * operators should be through the get() methods. You will need to write a
 * couple of create() methods, but you'll probably never need
 * do call one.
 *
 * The get methods don't construct operators directly; if the operator
 * already exists, a cached value is returned. Otherwise, the create()
 * method is called and the newly-created operator is cached. 
 * 
 * <H1> How to implement a Hessian factory for your problem </H1>
 *
 * There are two pure virtual methods you must implement: createA() 
 * and createC(), which produce the single-timestep advance operator
 * and observation operator, respectively. 
 * 
 * The other create() methods have default implementations in terms
 * of implicit operators, and normally won't need to be overridden. 
 * In some cases you may want to override them; for example, you might
 * want to use an explicitly-formed transpose or matrix-matrix product
 * rather than the implicit implementations of these. 
 */
template <class Scalar> 
class LTIProblemFactoryBase
{
public:
  /** Constructor */
  LTIProblemFactoryBase(int nSteps)
    : nSteps_(nSteps), 
      A_(), C_(), 
      At_(), Ct_(), 
      bigAInv_(), bigAInvT_(),
      bigF_(), bigFt_(),
      bigC_(), bigCt_(),
      H_()
    {}

  /** */
  virtual ~LTIProblemFactoryBase() {;}

  /** \name Access to single-timestep operators */
  //@{
  /** Access the single-timestep advance operator */
  LinearOperator<Scalar> getA() const 
    {
      if (A_.ptr().get()==0) A_ = createA();
      return A_;
    }

  /** Access the single-timestep observation operator */
  LinearOperator<Scalar> getC() const
    {
      if (C_.ptr().get()==0) C_ = createC();
      return C_;
    }


  /** Access the single-timestep adjoint advance operator */
  LinearOperator<Scalar> getAt() const 
    {
      if (At_.ptr().get()==0) At_ = createAt();
      return At_;
    }


  /** Access the single-timestep adjoint observation operator */
  LinearOperator<Scalar> getCt() const 
    {
      if (Ct_.ptr().get()==0) Ct_ = createCt();
      return Ct_;
    }
  //@}

  /** \name Multiple-timestep operators */
  //@{
  /** Access the multiple-timestep advance operator */
  LinearOperator<Scalar> getBigAInv() const 
    {
      if (bigAInv_.ptr().get()==0) bigAInv_ = createBigAInv();
      return bigAInv_;
    }

  /** Access the zero-padding operator */
  LinearOperator<Scalar> getBigF() const 
    {
      if (bigF_.ptr().get()==0) bigF_ = createBigF();
      return bigF_;
    }

  /** Access the multiple-timestep observation operator */
  LinearOperator<Scalar> getBigC() const
    {
      if (bigC_.ptr().get()==0) bigC_ = createBigC();
      return bigC_;
    }


  /** Access the multiple-timestep adjoint advance operator */
  LinearOperator<Scalar> getBigAInvT() const 
    {
      if (bigAInvT_.ptr().get()==0) bigAInvT_ = createBigAInvT();
      return bigAInvT_;
    }


  /** Access the multiple-timestep adjoint observation operator */
  LinearOperator<Scalar> getBigCt() const 
    {
      if (bigCt_.ptr().get()==0) bigCt_ = createBigCt();
      return bigCt_;
    }


  /** Access the multiple-timestep adjoint zero-padding operator */
  LinearOperator<Scalar> getBigFt() const 
    {
      if (bigFt_.ptr().get()==0) bigFt_ = createBigFt();
      return bigFt_;
    }

  /** Access the space-time Hessian. This consists of a forward
   * calculation followed by an adjoint calculation */
  LinearOperator<Scalar> getH() const 
    {
      if (H_.ptr().get()==0) H_ = createH();
      return H_;
    }
  //@}

protected:

  /** \name Pure virtual functions you must implement when setting up a concrete LTI problem */
  //@{
  /** Create the operator that advances the system through one timestep. */
  virtual LinearOperator<Scalar> createA() const = 0 ;

  /** Create the operator that produces the observable quantities given
   * a state. */
  virtual LinearOperator<Scalar> createC() const = 0 ;
  //@}

  /** \name Utilities */
  //@{
  /** Create a block vector space in which a space is replicated n times. */
  VectorSpace<Scalar> blockSpace(
    int n,
    const VectorSpace<Scalar>& sp
    ) const 
    {
      Array<VectorSpace<Scalar> > s(n, sp);
      return productSpace(s);
    }

  /** Create a block diagonal operator with a single block repeated
   * on all diagonal entries. */
  LinearOperator<Scalar> blockDiag(
    int n,
    const LinearOperator<Scalar>& C
    ) const 
    {
      Array<VectorSpace<Scalar> > d(n, C.domain());
      Array<VectorSpace<Scalar> > r(n, C.range());
      
      RCP<LinearOpBase<Scalar> > op
        = rcp(new BlockOperator<Scalar>(productSpace(d), productSpace(r)));
      LinearOperator<Scalar> rtn = op;
      for (int i=0; i<n; i++) rtn.setBlock(i, i, C);
      rtn.endBlockFill();
      return rtn;
    }
  //@}
  
  /** \name Create() methods with good default implementations */
  //@{
  /** */
  virtual LinearOperator<Scalar> createBigAInv() const 
    {
      RCP<LinearOpBase<Scalar> > rtn 
        = rcp(new InverseLTIOp<Scalar>(nSteps_, getA(), getAt()));
      return rtn;
    }

  /** */
  virtual LinearOperator<Scalar> createBigC() const 
    {return blockDiag(nSteps_, this->getC());}

  /** */
  virtual LinearOperator<Scalar> createBigF() const 
    {
      LinearOperator<Scalar> A = this->getA();
      std::cout << "A.domain().dim() = " << A.domain().dim() << std::endl;
      VectorSpace<Scalar> littleDomain = productSpace<Scalar>(tuple(A.domain()));
      
      LinearOperator<Scalar> I = identityOperator<Scalar>(A.domain());
      LinearOperator<Scalar> Z = zeroOperator<Scalar>(A.domain(), A.range());
      VectorSpace<Scalar> bigRange = this->blockSpace(nSteps_, A.range());
      std::cout << "bigRange.dim() = " << bigRange.dim() << std::endl;

      RCP<LinearOpBase<Scalar> > op
        = rcp(new BlockOperator<Scalar>(littleDomain, bigRange));
      LinearOperator<Scalar> rtn = op;

      rtn.setBlock(0, 0, I);
      for (int i=1; i<bigRange.numBlocks(); i++)
      {
        rtn.setBlock(i, 0, Z);
      }
      rtn.endBlockFill();
      return rtn;
    }

  /** By default, we create an implicit transpose of A. However, when
   * a transpose solver is unavailable, we should override this
   * function with one that returns an explicit transpose
   * of A. See the class ExplicitlyTransposedLTIProblemFactory for
   * a simple implementation of this. */
  virtual LinearOperator<Scalar> createAt() const 
    {return this->getA().transpose();}

  /** */
  virtual LinearOperator<Scalar> createCt() const 
    {return this->getC().transpose();}

  /** */
  virtual LinearOperator<Scalar> createBigAInvT() const 
    {return this->getBigAInv().transpose();}

  /** */
  virtual LinearOperator<Scalar> createBigFt() const 
    {return this->getBigF().transpose();}

  /** */
  virtual LinearOperator<Scalar> createBigCt() const 
    {return this->getBigC().transpose();}

  /** */
  virtual LinearOperator<Scalar> createH() const 
    {
      LinearOperator<Scalar> bigAInv = this->getBigAInv();
      LinearOperator<Scalar> bigAInvT = this->getBigAInvT();
      LinearOperator<Scalar> bigF = this->getBigF();
      LinearOperator<Scalar> bigC = this->getBigC();
      LinearOperator<Scalar> bigFt = this->getBigFt();
      LinearOperator<Scalar> bigCt = this->getBigCt();

      H_ = bigFt*bigAInvT*bigCt*bigC*bigAInv*bigF;
      return H_;
    }
  //@}
  
private:
  int nSteps_;
  mutable LinearOperator<Scalar> A_;
  mutable LinearOperator<Scalar> C_;

  mutable LinearOperator<Scalar> At_;
  mutable LinearOperator<Scalar> Ct_;

  mutable LinearOperator<Scalar> bigAInv_;
  mutable LinearOperator<Scalar> bigAInvT_;

  mutable LinearOperator<Scalar> bigF_;
  mutable LinearOperator<Scalar> bigFt_;

  mutable LinearOperator<Scalar> bigC_;
  mutable LinearOperator<Scalar> bigCt_;

  mutable LinearOperator<Scalar> H_;
};
}


#endif
