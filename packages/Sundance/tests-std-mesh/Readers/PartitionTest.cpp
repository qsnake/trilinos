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

class PartitionField : public FieldBase
{
public: 
  /** */
  PartitionField(bool isElem, int dim, const Array<int>& assignments)
    : isElem_(isElem), dim_(dim), data_(assignments)
    {}

  /** */
  virtual ~PartitionField() {;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const 
    {
      TEST_FOR_EXCEPT(!isDefined(cellDim, cellID, elem));
      return data_[cellID];
    }
  
  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const 
    {
      bool dimOK = (isElem_ && (cellDim == dim_)) || (!isElem_ && cellDim==0);
        
      return dimOK && cellID < (int) data_.size();
    }
  
  /** */
  virtual bool isPointData() const {return !isElem_;}

  /* */
  GET_RCP(FieldBase);

private:
  bool isElem_;
  int dim_;
  Array<int> data_;
  
};

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusMeshReader("../../../tests-std-mesh/Readers/wheel", meshType);

      Mesh mesh = mesher.getMesh();

      int np = 4;

      RCP<SerialPartitionerBase> part 
        = rcp(new FileIOChacoPartitioner("part"));

      Array<Mesh> submesh = part->makeMeshParts(mesh, np);


      for (int p=0; p<np; p++)
      {
        std::string filename = "wheel";
        FieldWriter w = new ExodusWriter(filename);
        w.impersonateParallelProc(np, p);
        w.addMesh(submesh[p]);
        w.write();
      }

      
      
      
      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << "Detected exception: " << e.what() << std::endl;
		}
}

