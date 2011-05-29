#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceMeshTransformation.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceExtrusionMeshTransformation.hpp"
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

      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 32, 1,
                                                         0.0, 1.0, 32, 1,
                                                         meshType);

      Mesh mesh2D = mesher.getMesh();

      MeshTransformation extruder = new ExtrusionMeshTransformation(0.0, 1.0, 32, meshType);

      Mesh mesh3D = extruder.apply(mesh2D);

      FieldWriter w3 = new VTKWriter("test3d");

      w3.addMesh(mesh3D);

      w3.write();

      std::cout << "num elements = " << mesh3D.numCells(3) << std::endl;
      std::cout << "num nodes = " << mesh3D.numCells(0) << std::endl;

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
}
