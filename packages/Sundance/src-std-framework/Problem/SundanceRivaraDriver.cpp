#include "SundanceRivaraDriver.hpp"
#include "SundanceMesh.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceExprFieldWrapper.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance::Rivara;
using Sundance::ExprFieldWrapper;


static Time& refTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mesh refinement"); 
  return *rtn;
}

static Time& m2rTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mesh to rivara copy"); 
  return *rtn;
}

static Time& r2mTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("rivara to mesh copy"); 
  return *rtn;
}

static Time& volSetTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("refinement stack building"); 
  return *rtn;
}


Mesh RefinementTransformation::apply(const Mesh& inputMesh) const 
{
  TimeMonitor timer(refTimer());

  int dim = inputMesh.spatialDim();
  MPIComm comm = inputMesh.comm();
  int numElems = inputMesh.numCells(dim);

  RCP<RivaraMesh> rivMesh = rcp(new RivaraMesh(dim, comm));

  Array<int> lidMap;

  meshToRivara(inputMesh, lidMap,rivMesh);
  
  ExprFieldWrapper expr(errExpr_);
  TEST_FOR_EXCEPTION(expr.isPointData(), RuntimeError,
    "Expected cell-based discrete function for area specification");

  {
    TimeMonitor timer(volSetTimer());
    numRefined_ = 0;
    for (int c=0; c<numElems; c++)
    {
      double err = expr.getData(dim,c,0);
      int rivLID = lidMap[c];
      Element* e = rivMesh->element(rivLID).get();
      double vol = e->volume();
      double newVol = vol * std::pow(reqErr_/(err+1.0e-12), 0.5*dim);
///      Out::os() << "c=" << c << " refine by " << newVol/vol << std::endl;
      if (newVol >= vol) continue;
      rivMesh->refinementSet().push(e);
      rivMesh->refinementAreas().push(newVol);
      numRefined_ ++;
    }
  }
  Out::os() << "refining n=" << numRefined_ << " cells" << std::endl;
  rivMesh->refine();


  Mesh outputMesh = rivaraToMesh(rivMesh, comm);

  return outputMesh;
}


void RefinementTransformation::meshToRivara(
  const Mesh& mesh, 
  Array<int>& lidMap,
  RCP<RivaraMesh>& rivMesh) const 
{
  TimeMonitor timer(m2rTimer());
  int dim = mesh.spatialDim();
  int numNodes = mesh.numCells(0);
  int numElems = mesh.numCells(dim);

  for (int n=0; n<numNodes; n++)
  {
    int gid = n;
    int label = mesh.label(0,n);
    int ownerProcID = mesh.ownerProcID(0,n);
    Point x = mesh.nodePosition(n);
    rivMesh->addVertex(gid, x, ownerProcID, label);
  }

  lidMap.resize(numElems);
  for (int e=0; e<numElems; e++)
  {
    int gid = e;
    int label = mesh.label(dim,e);
    int ownerProcID = mesh.ownerProcID(dim,e);
    Array<int> verts;
    Array<int> fo;
    mesh.getFacetArray(dim, e, 0, verts, fo);
    int elemLID = rivMesh->addElement(gid, verts, ownerProcID, label);
    lidMap[e] = elemLID;
    /* label edges or faces */
    if (dim==2)
    {
      Array<int> edgeLIDs;
      mesh.getFacetArray(dim, e, 1, edgeLIDs, fo);
      for (int f=0; f<3; f++)
      {
        int edgeLabel = mesh.label(1,edgeLIDs[f]);
        rivMesh->element(elemLID)->edge(f)->setLabel(edgeLabel);
      }
    }
    else if (dim==3)
    {
      Array<int> faceLIDs;
      mesh.getFacetArray(dim, e, 2, faceLIDs, fo);
      for (int f=0; f<4; f++)
      {
        int faceLabel = mesh.label(2,faceLIDs[f]);
        rivMesh->element(elemLID)->face(f)->setLabel(faceLabel);
      }
    }
  }
}


Mesh RefinementTransformation::rivaraToMesh(const RCP<RivaraMesh>& rivMesh,
  const MPIComm& comm) const 
{
  TimeMonitor timer(r2mTimer());
  int dim = rivMesh->spatialDim();
  int numNodes = rivMesh->numNodes();

  Mesh mesh = meshType_.createEmptyMesh(dim, comm);

  for (int n=0; n<numNodes; n++)
  {
    const RCP<Node>& node = rivMesh->node(n);
    const Point& x = node->pt();
    int gid = node->globalIndex();
    int ownerProcID = node->ownerProc();
    int label = node->label();
    mesh.addVertex(gid, x, ownerProcID, label);
  }


  ElementIterator iter(rivMesh.get());

  int gid=0;

  Array<int> verts(dim+1);
  Array<int> fo;
  Array<int> edgeLIDs;
  Array<int> faceLIDs;
      
  while (iter.hasMoreElements())
  {
    const Element* e = iter.getNextElement();
    int ownerProcID = e->ownerProc();
    int label = e->label();
    const Array<RCP<Node> >& nodes = e->nodes();

    for (int i=0; i<nodes.size(); i++)
    {
      verts[i] = nodes[i]->globalIndex();
    }
    int lid = mesh.addElement(gid, verts, ownerProcID, label);
    gid++;
    /* label edges or faces */
    if (dim==2)
    {
      if (e->hasNoEdgeLabels()) continue;
      mesh.getFacetArray(dim, lid, 1, edgeLIDs, fo);
      for (int f=0; f<3; f++)
      {
        int edgeLabel = e->edge(f)->label();
        mesh.setLabel(1, edgeLIDs[f], edgeLabel);
      }
    }
    else if (dim==3)
    {
      if (e->hasNoFaceLabels()) continue;
      mesh.getFacetArray(dim, lid, 2, faceLIDs, fo);
      for (int f=0; f<4; f++)
      {
        int faceLabel = e->face(f)->label();
        mesh.setLabel(2, faceLIDs[f], faceLabel);
      }
    }
  }

  return mesh;
  
}

  
