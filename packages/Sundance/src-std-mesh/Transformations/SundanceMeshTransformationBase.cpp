#include "SundanceMeshTransformationBase.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

Mesh MeshTransformationBase::createMesh(int dim, const MPIComm& comm) const 
{
  return meshType_.createEmptyMesh(dim, comm);
}


