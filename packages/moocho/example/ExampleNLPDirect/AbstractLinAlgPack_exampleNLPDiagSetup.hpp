// ////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_exampleNLPDiagSetup.hpp

#ifndef ALAP_EXPL_NLP_DIAG_SETUP_HPP
#define ALAP_EXPL_NLP_DIAG_SETUP_HPP

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief Create a vector space given the input arguments argc, argv[] and an MPI communicator
 *
 */
int exampleNLPDiagSetup(
  int argc, char* argv[], MPI_Comm comm
  ,Teuchos::RCP<const VectorSpace> *vec_space
  ,int *n, value_type *xo, bool *has_bounds, bool *dep_bounded
  );

} // namespace AbstractLinAlgPack

#endif // ALAP_EXPL_NLP_DIAG_SETUP_HPP
