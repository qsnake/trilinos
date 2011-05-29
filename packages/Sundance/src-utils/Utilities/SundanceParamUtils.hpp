/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_PARAMUTILS_H
#define SUNDANCE_PARAMUTILS_H

#include "SundanceDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Sundance
{
Teuchos::ParameterList mergeParamLists(
  const Teuchos::ParameterList& pDefault, 
  const Teuchos::ParameterList& pIn);
}

#endif
