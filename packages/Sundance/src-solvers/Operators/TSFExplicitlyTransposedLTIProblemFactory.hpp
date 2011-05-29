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

#ifndef TSF_EXPLICITLYTRANSPOSED_LTI_PROBLEMFACTORY_HPP
#define TSF_EXPLICITLYTRANSPOSED_LTI_PROBLEMFACTORY_HPP


#include "SundanceDefs.hpp"
#include "TSFDefaultLTIProblemFactory.hpp"
#include "TSFEpetraMatrix.hpp"
#include "EpetraExt_Transpose_RowMatrix.h"


namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;

/** 
 * This class produces LTI problems where any transpose solves are
 * done using explicit transforms. For example, with backwards Euler
 * time stepping, the advance matrix \f$ A\f$ will be 
 * \f$ A = M B^{-1}, \f$ 
 * where \f$ M \f$ is the mass matrix and \f$ B = M - \Delta t K \f$. However,
 * with the AztecOO solver we can't do a transpose solve on B, so we will want
 * to form \f$ B^T\f$ explicitly, after which we can form \f$ A^T \f$ 
 * implicitly as \f$ A^T = (B^T)^{-1} M^T.\f$ In this formulation
 * we are only doing a solve on the {\it explicit} transpose of \f$ B. \f$
 */
template <class Scalar> 
class ExplicitlyTransposedLTIProblemFactory
  : public DefaultLTIProblemFactory<Scalar>
{
public:
  /** Constructor */
  ExplicitlyTransposedLTIProblemFactory(int nSteps)
    : DefaultLTIProblemFactory<Scalar>(nSteps),
      inputAt_()
    {}



  /** 
   * After construction, the user should set the operators by calling
   * init().
   *
   * \param A The "small" A operator. Small A does the advance through 
   * a single timestep.
   * \param C The "small" C operator. Small C does state observations
   * at a single timestep.
   * \param At The explicit transpose of small A. 
   */
  void init(
    const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& At,
    const LinearOperator<Scalar>& C
    ) 
    {
      /* Call the init method on the default type. This will initialize
       * the "small" A and C operators. */
      DefaultLTIProblemFactory<Scalar>::init(A, C);

      /* Set the "small" A^T operator to the value given by the user */
      inputAt_ = At;
    }

protected:    

  /** This is a helper function that forms an explicit transpose of
   * an operator. For now, this only works if the operator is an 
   * Epetra_CrsMatrix. We'll want to generalize this later to other operator
   * types.
   *
   * One clumsy feature of this function is that because of the way the 
   * EpetraExt explicit transposer has been written, the storage for
   * the transposed operator is managed by the transposer object. 
   * Therefore, we need to keep a copy of the transposer object around 
   * for as long as the transposed matrix exists. We keep a copy
   * of it here, as the transposer_  data member of the factory object. 
   * This means that the factory object had better not go out of scope while
   * the matrices it produces still exist. This is dangerous. We should
   * talk with Rob H and Mike H about redesigning the EpetraExt transformation
   * objects so that the objects they produce have integrity independent
   * of the transformation objects.
   */
  LinearOperator<Scalar> formExplicitTranspose(const LinearOperator<Scalar>& X) const 
    {
      /* create a transposer object, and store it in the transposer_ 
       * data member. See the comment above about why ythis is bad but
       * necessary. */
      transposer_ = rcp(new EpetraExt::RowMatrix_Transpose(true));

      /* The transposer will act on an epetra crs matrix. Extract a 
       * Crs matrix from the input operator X. This will throw an exception
       * if the object handled by X is not an EpetraMatrix. */
      Epetra_CrsMatrix& epX = EpetraMatrix::getConcrete(X);

      /* Do the transpose operation, which returns a Crs matrix but in the
       * form of its base class, an Epetra_RowMatrix.  */
      Epetra_RowMatrix& eprXt = (*transposer_)(epX);

      /* The TSF EpetraMatrix works with Epetra_CrsMatrix, so we need
       * to cast the RowMatrix to a CrsMatrix. This should always work, 
       * because the transposer is implemented in terms of CrsMatrix. If it
       * doesn't work, there's been an error somehere in Trilinos. */
      Epetra_CrsMatrix* epXt = dynamic_cast<Epetra_CrsMatrix*>(&eprXt);
      TEST_FOR_EXCEPTION(epXt == 0, std::runtime_error, "expected return type "
        "of EpetraExt tranposer to be a CrsMatrix");

      
      /* We need Epetra spaces for the domain and range spaces of the new
       * operator. We simply get these from the original, untransposed 
       * operator and then swap the order of domain and range. The operator
       * interface returns generic Thyra VectorSpace objects, which we
       * dynamic cast into Epetra VectorSpace objects. This cast should
       * always work. 
       */ 
      RCP<const EpetraVectorSpace> epXRange = rcp_dynamic_cast<const EpetraVectorSpace>(X.range().ptr());
      RCP<const EpetraVectorSpace> epXDomain = rcp_dynamic_cast<const EpetraVectorSpace>(X.domain().ptr());


      /* We can now create a TSF linear operator for the transpose. We need 
       * to pass an ownership flag of "false" for the RCP of the 
       * Epetra_CrsMatrix because the transposed matrix is
       * owned by the transposer. We also swap the original operator's 
       * domain and range spaces in creating the transpose. */
      RCP<LinearOpBase<Scalar> > XtPtr = rcp(new EpetraMatrix(rcp(epXt, false), 
          epXRange, epXDomain));
      LinearOperator<Scalar> Xt = XtPtr;

      /** all done! */
      return Xt;
    }

  /** Return the transpose of A */
  virtual LinearOperator<Scalar> createAt() const 
    {
      TEST_FOR_EXCEPT(inputAt_.ptr().get() == 0);
      return inputAt_;
    }

private:
  LinearOperator<Scalar> inputAt_;
  /** Hack hack hack.... we need to keep the Epetra_Transposer object around
   * because it manages the memory of the transposed matrix. */
  mutable RCP<EpetraExt::RowMatrix_Transpose> transposer_;
};
}


#endif
