// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability 
//                   of Abstract Numerical Algorithms 
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
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_RTOP_SUNDIALS_HELPERS_HPP
#define RTOPPACK_RTOP_SUNDIALS_HELPERS_HPP

#include "RTOpPack_RTOpTHelpers.hpp"



namespace RTOpPack 
{
  /**
   *
   */
#ifdef TRILINOS_6.0

  template<class Scalar>
  class RTOpBoolReduceAndTransform 
    : public ROpIndexReductionBase<Scalar>,
      public ROpScalarTransformationBase<Scalar>
{

  public:
    typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
    /** */
    RTOpBoolReduceAndTransform()
      : RTOpT<Scalar>(""), 
        ROpIndexReductionBase<Scalar>(1),
        ROpScalarTransformationBase<Scalar>() 
    {;}

    /** */
    virtual ~RTOpBoolReduceAndTransform(){;}

    /** */
    index_type operator()(const ReductTarget& reduct_obj ) const 
    { return this->getRawVal(reduct_obj); }


    /** \brief Default implementation here is for a logical AND. */
    void reduce_reduct_objs(const ReductTarget& in_reduct_obj, 
                            ReductTarget* inout_reduct_obj) const
    {
      const index_type in_val    = this->getRawVal(in_reduct_obj);
      const index_type inout_val = this->getRawVal(*inout_reduct_obj);
      this->setRawVal( in_val && inout_val, inout_reduct_obj );
    }

  };

#endif /* TRILINOS_6 */

}

/** \brief Use within an apply_op(...) function implementation where num_vecs==3, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_3_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION(                                                   \
                     (NUM_VECS)!=3 || (SUB_VECS)==NULL                  \
                     ,RTOpPack::InvalidNumVecs                          \
                     ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==3, sub_vecs!=NULL" \
                     );                                                 \
  TEST_FOR_EXCEPTION(                                                   \
                     (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL        \
                     ,RTOpPack::InvalidNumTargVecs                      \
                     ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
                     );                                                 \
  TEST_FOR_EXCEPTION(                                                   \
                     (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim()   \
                     || (SUB_VECS)[0].subDim() != (SUB_VECS)[2].subDim() \
                     ||(SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
                     ||(SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
                     ,IncompatibleVecs                                  \
                     ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
                     );                                                 \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim();        \
  const Scalar                 *v0_val = (SUB_VECS)[0].values();        \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride();        \
  const Scalar                 *v1_val = (SUB_VECS)[1].values();        \
  const ptrdiff_t              v1_s    = (SUB_VECS)[1].stride();        \
  const Scalar                 *v2_val = (SUB_VECS)[2].values();        \
  const ptrdiff_t              v2_s    = (SUB_VECS)[2].stride();


#endif
