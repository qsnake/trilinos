#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"


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



int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);



      Mesh mesh = mesher.getMesh();

      FieldWriter w = new VerboseFieldWriter();

      w.addMesh(mesh);

      w.write();

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
}
