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

#ifndef TSFINVERSELTIOP_HPP
#define TSFINVERSELTIOP_HPP

#include "SundanceDefs.hpp"
#include "TSFSimplifiedLinearOpBase.hpp"
#include "TSFHomogeneouslyBlockedLinearOp.hpp"
#include "TSFSimpleIdentityOpDecl.hpp"


namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

/** 
 * InverseLTIOp encapsulates a block operator representation
 * of timestepping for a linear, time-invariant system of ODEs.
 * The ODEs are replaced by a discrete-time system
 * \f[
 * {\bf x}_{k+1} = {\bf A} \cdot {\bf x}_k.
 * \f]
 * This can be written as the block system
 * \f[
 * \left[
 * \begin{array}{cccc}
 * I  & 0  & 0  & \cdots \\
 * -A & I  & 0  & \\
 * 0  & -A & I  & \\
 * \vdots & & & \ddots 
 * \end{array}
 * \right]
 * \left[\begin{array}{c} x_0 \\ x_1 \\ x_2 \\ \vdots \end{array}\right]
 * =
 * \left[\begin{array}{c} x_0 \\ 0 \\ 0 \\ \vdots \end{array}\right]
 * \f].
 * The inner matrix \f$ A\f$ depends on the timestepping algorithm
 * used as well as the system of ODEs. 
 *
 * This class solves the system by backsubstitution.
 */
template <class Scalar> 
class InverseLTIOp
  : public virtual LinearOpBase<Scalar>,
    public virtual HomogeneouslyBlockedLinearOp<Scalar>,
    public virtual SimplifiedLinearOpBase<Scalar>
{
public:

  /** 
   * Construct a InverseLTIOp that takes <t>numTimesteps</t> steps
   * with the operator \f$ A\f$.
   */
  InverseLTIOp(int numTimesteps, const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& At)
    : HomogeneouslyBlockedLinearOp<Scalar>(
      A.domain(), numTimesteps,
      A.range(), numTimesteps),
      A_(A), At_(At)
    {}



  /** 
   * Apply the operator
   */
  void applyOp(
    const Thyra::EOpTransp M_trans,
    const Vector<Scalar>& in,
    Vector<Scalar> out
    ) const 
    {
      int nbIn = in.space().numBlocks();
      int nbOut = out.space().numBlocks();
      TEST_FOR_EXCEPTION(nbIn != nbOut, std::runtime_error,
        "expected a square block structure, found nbIn=" << nbIn
        << ", nbOut=" << nbOut);
          

      LinearOperator<Scalar> I = identityOperator<Scalar>(in.space().getBlock(0));

      if (M_trans==Thyra::NOTRANS)
      {
        for (int i=0; i<this->numBlockRows(); i++)
        {
          if (i==0)
          {
            out.setBlock(i, I*in.getBlock(i));
          }
          else
          {
            Vector<Scalar> xi1 = out.getBlock(i-1);
            Vector<Scalar> bi = in.getBlock(i);
            out.setBlock(i, bi + A_*xi1);
          }
        }
      }
      else if (M_trans==Thyra::TRANS)
      {
        for (int i=this->numBlockCols()-1; i>=0; i--)
        {
          if (i==this->numBlockCols()-1)
          {
            out.setBlock(i, I*in.getBlock(i));
          }
          else
          {
            Vector<Scalar> bi = in.getBlock(i);
            Vector<Scalar> xi1 = out.getBlock(i+1);
            out.setBlock(i, bi + At_*xi1);
          }
        }
      }
      else
      {
        TEST_FOR_EXCEPT(true);
      }

    }

private:
  LinearOperator<Scalar> A_;
  LinearOperator<Scalar> At_;

};
}


#ifdef TRILINOS_6
#undef DefaultColumnwiseMultiVector
#endif

#endif
