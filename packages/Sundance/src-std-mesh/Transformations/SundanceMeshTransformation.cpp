#include "SundanceMeshTransformation.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;



MeshTransformation::MeshTransformation()
  : Handle<MeshTransformationBase>()
{}

MeshTransformation::MeshTransformation(Handleable<MeshTransformationBase>* rawPtr)
  : Handle<MeshTransformationBase>(rawPtr)
{}


MeshTransformation::MeshTransformation(const RCP<MeshTransformationBase>& smartPtr)
  : Handle<MeshTransformationBase>(smartPtr)
{}

Mesh MeshTransformation::apply(const Mesh& inputMesh) const 
{
  Mesh rtn = ptr()->apply(inputMesh);
  //if (rtn.spatialDim() > 1) rtn.assignIntermediateCellOwners(1);
  //if (rtn.spatialDim() > 2) rtn.assignIntermediateCellOwners(2);
  return rtn;
}
