#include "SundanceMeshSourceBase.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


MeshSourceBase::MeshSourceBase(const MeshType& meshType,
                               const MPIComm& comm)
  : cachedMesh_(),
    hasCachedMesh_(),
    meshType_(meshType),
    comm_(comm),
    elemAttributes_(rcp(new Array<Array<double> >())),
    nodeAttributes_(rcp(new Array<Array<double> >()))
{
}

MeshSourceBase::MeshSourceBase(const ParameterList& params)
  : cachedMesh_(),
    hasCachedMesh_(),
    meshType_(new BasicSimplicialMeshType()),
    comm_(MPIComm::world()),
    elemAttributes_(rcp(new Array<Array<double> >())),
    nodeAttributes_(rcp(new Array<Array<double> >()))
{
  
}

Mesh MeshSourceBase::getMesh() const
{
  Tabs tabs;
  
  /* if we don't have a cached mesh, build one */
  if (!hasCachedMesh_)
    {
      Mesh rtn =  fillMesh();
      if (verb() > 0)
        {
          std::cerr << tabs << "got a mesh with " << rtn.numCells(0)
               << " nodes and " << rtn.numCells(rtn.spatialDim())
               << " maximal cells" << std::endl;
        }
      return rtn;
    }
  return cachedMesh_;
}

void MeshSourceBase
::getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
                RCP<Array<Array<double> > >& elemAttributes) const
{
  nodeAttributes = nodeAttributes_;
  elemAttributes = elemAttributes_;
}

Mesh MeshSourceBase::createMesh(int dim) const
{
  cachedMesh_ = meshType_.createEmptyMesh(dim, comm_);
  cachedMesh_.ptr()->setVerbosity(verb());
  hasCachedMesh_ = true;
  
  return cachedMesh_;
}
