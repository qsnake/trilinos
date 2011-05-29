// ////////////////////////////////////////////////////
// AbstractLinAlgPack_ThyraAccessors.hpp

#ifndef ALAP_THYRA_ACCESSORS_HPP
#define ALAP_THYRA_ACCESSORS_HPP

#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"

namespace AbstractLinAlgPack {

/** \brief . */
void get_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,const Vector                                                  &vec
  ,Teuchos::RCP<const Thyra::VectorBase<value_type> >    *thyra_vec
  );

/** \brief . */
void free_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,const Vector                                                  &vec
  ,Teuchos::RCP<const Thyra::VectorBase<value_type> >    *thyra_vec
  );

/** \brief . */
void get_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,VectorMutable                                                 *vec
  ,Teuchos::RCP<Thyra::VectorBase<value_type> >          *thyra_vec
  );

/** \brief . */
void commit_thyra_vector(
  const VectorSpaceThyra                                         &thyra_vec_spc
  ,VectorMutable                                                 *vec
  ,Teuchos::RCP<Thyra::VectorBase<value_type> >          *thyra_vec
  );

}

#endif // ALAP_THYRA_ACCESSORS_HPP
