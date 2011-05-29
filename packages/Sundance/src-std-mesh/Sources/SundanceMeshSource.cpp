#include "SundanceMeshSource.hpp"
#include "SundanceOut.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


static Time& getMeshTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("get mesh"); 
  return *rtn;
}

MeshSource::MeshSource()
  : Handle<MeshSourceBase>()
{}

MeshSource::MeshSource(Handleable<MeshSourceBase>* rawPtr)
  : Handle<MeshSourceBase>(rawPtr)
{}


MeshSource::MeshSource(const RCP<MeshSourceBase>& smartPtr)
  : Handle<MeshSourceBase>(smartPtr)
{}

Mesh MeshSource::getMesh() const
{
  TimeMonitor timer(getMeshTimer());

  Mesh rtn;
  try
    {
      Tabs tabs;
      int nProc = ptr()->comm().getNProc();
      SUNDANCE_OUT(ptr()->verb() > 0, 
                   "MeshSource::getMesh()");
      if (staggerOutput() && nProc > 1)
        {
          int myRank = ptr()->comm().getRank();
          for (int p=0; p<nProc; p++)
            {
              ptr()->comm().synchronize();
              if (p != myRank) continue;
              SUNDANCE_OUT(ptr()->verb() > 0, 
                           "========= Building local mesh on processor " 
                           << p << " ============ ");
              rtn = ptr()->getMesh();
            }
        }
      else 
        {
          rtn = ptr()->getMesh();
        }

      if (rtn.spatialDim() > 1) rtn.assignIntermediateCellGIDs(1);
      if (rtn.spatialDim() > 2) rtn.assignIntermediateCellGIDs(2);
    }
  catch(std::exception& e)
    {
      SUNDANCE_TRACE(e);
    }
  return rtn;
}

void MeshSource::getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
                               RCP<Array<Array<double> > >& elemAttributes) const
{
  getMesh();
  ptr()->getAttributes(nodeAttributes, elemAttributes);
}


MeshType& MeshSource::defaultMeshType() 
{
  static MeshType rtn = new BasicSimplicialMeshType(); 
  return rtn;
}

const MPIComm& MeshSource::comm() const
{
  return ptr()->comm();
}
