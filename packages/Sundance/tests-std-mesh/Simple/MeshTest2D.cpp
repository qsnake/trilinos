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



int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      int np = MPIComm::world().getNProc();

      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 4, 1,
                                                         0.0, 1.0, 4, np,
                                                         meshType);



      Mesh mesh = mesher.getMesh();

      FieldWriter w = new VerboseFieldWriter();
      FieldWriter w1 = new VerboseFieldWriter("test2d");
      FieldWriter w2 = new TriangleWriter("test2d");
      FieldWriter w3 = new VTKWriter("test2d");

      mesh.verbosity() = VerbExtreme;

      w.addMesh(mesh);
      w1.addMesh(mesh);
      w2.addMesh(mesh);
      w3.addMesh(mesh);

      w.write();
      w1.write();
      w2.write();
      w3.write();
      TimeMonitor::summarize();
      
    }
	catch(std::exception& e)
		{
      std::cerr << "Detected exception: " << e.what() << std::endl;
		}
}
