// /////////////////////////////////////////////////////
// VectorApplyOpSerialBase.cpp

#include "AbstractLinAlgPack_VectorApplyOpSerialBase.hpp"
#include "AbstractLinAlgPack_apply_op_helper.hpp"

namespace AbstractLinAlgPack {

VectorApplyOpSerialBase::VectorApplyOpSerialBase()
  : in_apply_op_(false)
{}

void VectorApplyOpSerialBase::apply_op_serial(
  const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
  ) const
{
  TEST_FOR_EXCEPTION(
    in_apply_op_, std::logic_error
    ,"VectorApplyOpSerialBase::apply_op_serial(...): Error, this function has been entered "
    "recursively which most likely means that the explicit sub-vector access methods Vector::get_sub_vector(...), "
    "Vector::free_sub_vector(...), VectorMutable::get_sub_vector(...), VectorMutable::commit_sub_vector(...) "
    "have not been overridden correctly on this concrete class \'" << typeName(*this) << "\' to not call "
    "apply_op(...) in there implemenations."
    );
  in_apply_op_ = true;
  AbstractLinAlgPack::apply_op_serial(
    op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
    ,first_ele,sub_dim,global_offset
    );
  in_apply_op_ = false;
}

} // namespace AbstractLinAlgPack
