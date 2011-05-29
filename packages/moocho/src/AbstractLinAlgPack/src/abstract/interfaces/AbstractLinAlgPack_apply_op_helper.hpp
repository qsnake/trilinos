// //////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_apply_op_helper.hpp

#ifndef APPLY_OP_HELPER_H
#define APPLY_OP_HELPER_H

#include "AbstractLinAlgPack_Types.hpp"
#include "RTOpPack_RTOpT.hpp"

namespace AbstractLinAlgPack {

/** \brief Validate the inputs to <tt>apply_op()</tt>.
 *
 * Throws an exception if one of the preconditions is not met.
 *
 * ToDo: Finish documentation.
 */
void apply_op_validate_input(
  const char func_name[]
  ,const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
  );

/** \brief Implements reduction/transformation operators for any serial vectors
 * using just the public vector interface.
 *
 * Note that this function does not validate the input arguments so it is up to
 * the client to do that (i.e. by calling \c apply_op_validate_input()).
 *
 * ToDo: Finish documentation!
 */
void apply_op_serial(
  const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
  );

} // end namespace AbstractLinAlgPack

#endif // APPLY_OP_HELPER_H
