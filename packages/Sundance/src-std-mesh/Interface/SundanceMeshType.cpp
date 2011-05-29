#include "SundanceMeshType.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

using std::endl;

MeshType::MeshType()
  : Handle<MeshTypeBase>()
{}

MeshType::MeshType(Handleable<MeshTypeBase>* rawPtr)
  : Handle<MeshTypeBase>(rawPtr)
{}


MeshType::MeshType(const RCP<MeshTypeBase>& smartPtr)
  : Handle<MeshTypeBase>(smartPtr)
{}

Mesh MeshType::createEmptyMesh(int dim, const MPIComm& comm) const 
{
  Mesh rtn;
  try
    {
      rtn = ptr()->createEmptyMesh(dim, comm);
    }
  catch(std::exception& e)
    {
      SUNDANCE_TRACE(e);
    }
  return rtn;
}


