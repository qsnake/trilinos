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
#include "SundanceFileIOChacoPartitioner.hpp"
#include "SundanceVTKWriter.hpp"
#include "Teuchos_CommandLineProcessor.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

using namespace std;

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

    CommandLineProcessor clp;
    TimeMonitor t(totalTimer());

    std::string infile="wheel";
    std::string outfile=infile;
    int numProc = 4;
    bool help=false;
    clp.setOption("i", &infile, "Input mesh filename");
    clp.setOption("o", &outfile, "Output mesh filename");
    clp.setOption("np", &numProc, "Number of partitions");
    clp.setOption("h", "nohelp", &help, "Help");
      
    clp.throwExceptions(false);

    CommandLineProcessor::EParseCommandLineReturn rtn 
      = clp.parse(argc, (char**) argv);

    TEST_FOR_EXCEPTION(rtn != CommandLineProcessor::PARSE_SUCCESSFUL,
      RuntimeError,
      "Command-line parsing failed");

    if (help)
    {
      cout << "Usage: partitionExo --i=inputFilename --o=outputFilename "
        "--np=numberOfPartitions" << std::endl;
      cout << "Do not include .exo suffix on filenames" << std::endl;
    }
    else
    {
      TEST_FOR_EXCEPT(infile.length()==0);
      TEST_FOR_EXCEPT(outfile.length()==0);
      TEST_FOR_EXCEPT(numProc<=1);


      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusMeshReader(infile, meshType);

      Mesh mesh = mesher.getMesh();

      RCP<SerialPartitionerBase> part 
        = rcp(new FileIOChacoPartitioner("part"));

      Array<Mesh> submesh = part->makeMeshParts(mesh, numProc);


      for (int p=0; p<numProc; p++)
      {
        FieldWriter w = new ExodusWriter(outfile);
        w.impersonateParallelProc(numProc, p);
        w.addMesh(submesh[p]);
        w.write();
      }
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    std::cerr << "Detected exception: " << e.what() << std::endl;
  }
}

