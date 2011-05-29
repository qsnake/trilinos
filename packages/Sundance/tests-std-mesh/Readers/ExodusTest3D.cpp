#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExodusMeshReader.hpp"
#include "SundanceExodusWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceVTKWriter.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;



static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


class ProcField : public FieldBase
{
public: 
  /** */
  ProcField(int dim)
    : dim_(dim), pid_(MPIComm::world().getRank())
    {}

  /** */
  virtual ~ProcField() {;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const 
    {
      TEST_FOR_EXCEPT(!isDefined(cellDim, cellID, elem));
      return pid_;
    }
  
  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const 
    {
      return cellDim == dim_;
    }
  
  /** */
  virtual bool isPointData() const {return dim_ == 0 ;}

  /* */
  GET_RCP(FieldBase);

private:
  int dim_;
  int pid_;
};

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusMeshReader("plate3D-0", meshType);



      Mesh mesh = mesher.getMesh();

      FieldWriter w = new VTKWriter("wheel");
      w.addMesh(mesh);
      w.addField("cell proc", new ProcField(mesh.spatialDim()));
      w.addField("point proc", new ProcField(0));
      w.write();

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << "Detected exception: " << e.what() << std::endl;
		}
}

