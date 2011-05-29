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

#ifndef TSFMULTIVECTOR_HPP
#define TSFMULTIVECTOR_HPP

#include "SundanceDefs.hpp"
#include "SundanceHandle.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFAccessibleVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Range1D.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFLinearOperatorDecl.hpp"



namespace TSFExtended
{
  using Thyra::Index;
  using Thyra::Range1D;
  using Thyra::EOpTransp;
  using namespace Teuchos;
  

  template <class Scalar>
  class MultiVector : public Sundance::Handle<Thyra::MultiVectorBase<Scalar> >
    {
    public:

      HANDLE_CTORS(MultiVector<Scalar>, Thyra::MultiVectorBase<Scalar>);
      
      /** Provide access to the columns as <\tt> Vector<\tt> objects */


      /** Get mutable view of column vector  */
      //      virtual Teuchos::RCP<Vector<Scalar> > col(Index j) const;
      virtual Vector<Scalar> col(Index j) const;



      /** Get nonmutable view of column vector  */
      //virtual Teuchos::RCP<const Vector<Scalar> > col(Index j) const;
      virtual const Vector<Scalar> col(Index j) const;



      /** Clone the multivector, if supported  */
      //virtual Teuchos::RCP<MultiVector<Scalar> > clone_mv() const;
      virtual MultiVector<Scalar> clone_mv() const;


        /** Returns a const sub-view of a contiguous set of columns of
	    the this multi-vector. */
//        virtual Teuchos::RCP<const MultiVector<Scalar> > 
//               subView( const Range1D& col_rng ) const;
       virtual const MultiVector<Scalar> subView(const Range1D& col_rng) const;


       /** Returns a non-const sub-view of a contiguous set of columns
	   of the t his multi-vector. */
//        virtual Teuchos::RCP<MultiVector<Scalar> > 
//               subView( const Range1D& col_rng );
       virtual MultiVector<Scalar> subView(const Range1D& col_rng);

        /** Returns a const sub-view of a non-contiguous set of columns of this
	    multi-vector. */
//       virtual Teuchos::RCP<const MultiVector<Scalar> > 
//             subView( const int numCols, const int cols[] ) const;
       virtual const MultiVector<Scalar> subView( const int numCols, 

						  const int cols[] ) const;
     

       /** Apply a reduction/transformation operator column by column
	   and return an array of the reduction objects.  */

//        virtual void applyOp(
//                 const RTOpPack::RTOpT<Scalar>   &primary_op
//                 ,const int                   num_multi_vecs
//                 ,const MultiVector<Scalar>*     multi_vecs[]
//                 ,const int                   num_targ_multi_vecs
//                 ,MultiVector<Scalar>*           targ_multi_vecs[]
//                 ,RTOpPack::ReductTarget*        reduct_objs[]
//                 ,const Index                    primary_first_ele
//                ,const Index                    primary_sub_dim
//                 ,const Index                    primary_global_offset
//                 ,const Index                    secondary_first_ele
//                 ,const Index                    secondary_sub_dim
//                 ) const;

       virtual void applyOp(
                const RTOpPack::RTOpT<Scalar>   &primary_op
                ,const int                   num_multi_vecs
                ,const MultiVector<Scalar>*     multi_vecs[]
                ,const int                   num_targ_multi_vecs
                ,MultiVector<Scalar>*           targ_multi_vecs[]
                ,RTOpPack::ReductTarget*        reduct_objs[]
                ,const Index                    primary_first_ele
               ,const Index                    primary_sub_dim
                ,const Index                    primary_global_offset
                ,const Index                    secondary_first_ele
                ,const Index                    secondary_sub_dim
                ) const;


        /** Apply a reduction/transformation operator column by column
            and reduce the intermediate reduction objects into one
            reduction object. */


       virtual void applyOp(
                const RTOpPack::RTOpT<Scalar>   &primary_op
                ,const RTOpPack::RTOpT<Scalar>  &secondary_op
                ,const int                   num_multi_vecs
                ,const MultiVector<Scalar>*     multi_vecs[]
                ,const int                   num_targ_multi_vecs
                ,MultiVector<Scalar>*           targ_multi_vecs[]
                ,RTOpPack::ReductTarget         *reduct_obj
                ,const Index                    primary_first_ele
                ,const Index                    primary_sub_dim
                ,const Index                    primary_global_offset
                ,const Index                    secondary_first_ele
                ,const Index                    secondary_sub_dim
                ) const;

      /** Get a non-mutable explicit view of a sub-multi-vector. */
         
       virtual void getSubMultiVector(
                const Range1D                       &rowRng
                ,const Range1D                      &colRng
                ,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
                ) const;



      /** Free an explicit view of a sub-multi-vector. */
       virtual void freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* 
					sub_mv ) const;




      /** Get a mutable explicit view of a sub-multi-vector. */
        virtual void getSubMultiVector(
                const Range1D                                &rowRng
                ,const Range1D                               &colRng
                ,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
                );




      /** Commit changes for a mutable explicit view of a sub-multi-vector. */
        virtual void commitSubMultiVector(
          RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv );


       /** @name Overridden methods from LinearOp */
        //@{

        ///
        /** This method is implemented in terms of the mulit-vector
	    <tt>applyOp( )</tt> method. */
       void apply( 
                const EOpTransp            M_trans 
                ,const Vector<Scalar>    &x 
                ,Vector<Scalar>          *y 
                ,const Scalar            alpha 
                ,const Scalar            beta 
                ) const; 

    /// This method is simply overridden to return <tt>this->clone_lons()</tt>.

       //Teuchos::RCP<const LinearOp<Scalar> > clone() const;
        LinearOperator<Scalar> clone() const;

        //@}
      
  };
  
}

 
#endif   
