// /////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorApplyOpSerialBase.hpp

#ifndef ALAP_VECTOR_APPLY_OP_SERIAL_BASE_HPP
#define ALAP_VECTOR_APPLY_OP_SERIAL_BASE_HPP

#include "AbstractLinAlgPack_Types.hpp"
#include "RTOpPack_RTOpT.hpp"

namespace AbstractLinAlgPack {

/** \brief This is a base class that is meant to be inherited by <tt>Vector</tt>
 * subclasses that wish to call <tt>TSFCore::apply_op_serial()</tt> when vectors
 * are in core.
 *
 * Calling this classes <tt>apply_op_serial()</tt> makes sure that the explicit
 * vector access functions have been implemented properly.
 */
class VectorApplyOpSerialBase {
public:

  /** \brief . */
  VectorApplyOpSerialBase();

  /** \brief . */
  void apply_op_serial(
    const RTOpPack::RTOp& op
    ,const size_t num_vecs, const Vector* vecs[]
    ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
    ,RTOpPack::ReductTarget *reduct_obj
    ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
    ) const;

private:

  mutable bool in_apply_op_;

}; // class VectorApplyOpSerialBase

} // namespace AbstractLinAlgPack

#endif // ALAP_VECTOR_APPLY_OP_SERIAL_BASE_HPP
