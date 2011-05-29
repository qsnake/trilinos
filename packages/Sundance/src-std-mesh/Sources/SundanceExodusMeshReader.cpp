#include "SundanceExodusMeshReader.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"
#include "SundancePathUtils.hpp"

#ifdef HAVE_SUNDANCE_EXODUS
#include "exodusII.h"
#endif 

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;
using namespace std;

ExodusMeshReader::ExodusMeshReader(const std::string& fname,
  const MeshType& meshType,
  const MPIComm& comm)
  : MeshReaderBase(fname, meshType, comm), 
    exoFilename_(fname),
    parFilename_(fname)
{
  if (nProc() > 1)
    {
      std::string suffix =  "-" + Teuchos::toString(nProc()) 
        + "-" + Teuchos::toString(myRank());
      exoFilename_ = exoFilename_ + suffix;
      parFilename_ = parFilename_ + suffix;
    }
  exoFilename_ = exoFilename_ + ".exo";
  parFilename_ = parFilename_ + ".pxo";
  
  setVerbosity( classVerbosity() );

  SUNDANCE_OUT(this->verb() > 1,
               "exodus filename = " << exoFilename_);
  
  if (nProc() > 1)
  {
    SUNDANCE_OUT(this->verb() > 1,
      "par filename = " << parFilename_);
  }
}





Mesh ExodusMeshReader::fillMesh() const 
{
  Mesh mesh;
#ifndef HAVE_SUNDANCE_EXODUS
  TEST_FOR_EXCEPTION(true, RuntimeError, 
    "ExodusMeshReader called for a build without ExodusII");
#else

  int CPU_word_size = 8;
  int IO_word_size = 0;
  float version;

  Array<int> ptGID;
  Array<int> ptOwner;
  Array<int> elemGID;
  Array<int> elemOwner;

  readParallelInfo(ptGID, ptOwner, elemGID, elemOwner);

  if (verb() > 2) ex_opts(EX_DEBUG | EX_VERBOSE);

  std::string resolvedName = searchForFile(exoFilename_);
  int exoID = ex_open(resolvedName.c_str(), EX_READ, 
    &CPU_word_size, &IO_word_size, &version);

  TEST_FOR_EXCEPTION(exoID < 0, RuntimeError, "ExodusMeshReader unable to "
    "open file: " << exoFilename_);

  TEST_FOR_EXCEPT(IO_word_size != 8 || CPU_word_size != 8);

  char title[MAX_LINE_LENGTH+1];

  int dim = 0;
  int numNodes = 0;
  int numElems = 0;
  int numElemBlocks = 0;
  int numNodeSets = 0;
  int numSideSets = 0;

  int ierr = ex_get_init(exoID, title, &dim, &numNodes, &numElems,
    &numElemBlocks, &numNodeSets, &numSideSets);

  TEST_FOR_EXCEPTION(numNodes <= 0, RuntimeError, "invalid numNodes=" 
    << numNodes);
  TEST_FOR_EXCEPTION(numElems <= 0, RuntimeError, "invalid numElems=" 
    << numElems);

  /* */
  if (nProc()==1)
  {
    ptGID.resize(numNodes);
    ptOwner.resize(numNodes);
    for (int i=0; i<numNodes; i++)
    {
      ptGID[i] = i;
      ptOwner[i] = 0;
    }
  }
  else
  {
    /* If we're running in parallel, we'd better have consistent numbers
     * of points in the .exo and .pex files. */
    TEST_FOR_EXCEPTION((int)ptGID.size() != numNodes, RuntimeError,
      "ExodusMeshReader::getMesh() found inconsistent "
      "numbers of points in .exo file and .pex files. Exo "
      "file " << exoFilename_ << " had nPoints=" 
      << numNodes << " but .pex file " 
      << parFilename_ << " had nPoints=" << ptGID.size());
  }

  /* now we can build the mesh */
  mesh = createMesh(dim);


  /* Read the points */
  Array<double> x(numNodes);
  Array<double> y(numNodes);
  Array<double> z(numNodes * (dim > 2));

  if (dim == 2)
  {
    ierr = ex_get_coord(exoID, &(x[0]), &(y[0]), (void*) 0);
    TEST_FOR_EXCEPT(ierr < 0);
    TEST_FOR_EXCEPT(ierr > 0);
  }
  else if (dim==3)
  {
    ierr = ex_get_coord(exoID, (void*) &(x[0]), (void*) &(y[0]), (void*) &(z[0]));
    TEST_FOR_EXCEPT(ierr < 0);
  }
  else 
  {
    TEST_FOR_EXCEPTION(dim < 2 || dim > 3, RuntimeError, 
      "invalid dimension=" << dim << " in ExodusMeshReader");
  }

  /* add the points to the mesh */
  for (int n=0; n<numNodes; n++)
  {
    Point p;
    if (dim==2)
    {
      p = Point(x[n], y[n]);
    }
    else
    {
      p = Point(x[n], y[n], z[n]);
    }
    mesh.addVertex(ptGID[n], p, ptOwner[n], 0);
  }


  /* Set up the element numbering */
  if (nProc()==1)
  {
    elemGID.resize(numElems);
    elemOwner.resize(numElems);
    for (int i=0; i<numElems; i++)
    {
      elemGID[i] = i;
      elemOwner[i] = 0;
    }
  }
  else
  {
    /* If we're running in parallel, we'd better have consistent numbers
     * of elements in the .exo and .pex files. */
    TEST_FOR_EXCEPTION((int)elemGID.size() != numElems, RuntimeError,
      "ExodusMeshReader::readElems() found inconsistent "
      "numbers of elements in .exo file and .pex files. Exodus "
      "file " << exoFilename_ << " had nElems=" 
      << numElems << " but .pex file " 
      << parFilename_ << " had nElems=" << elemGID.size());
  }

  /* Read the elements for each block */

  Array<int> blockIDs(numElemBlocks);
  if (numElemBlocks > 0)
  {
    ierr = ex_get_elem_blk_ids(exoID, &(blockIDs[0]));
  }
  TEST_FOR_EXCEPT(ierr < 0);
  int count = 0;
  for (int b=0; b<numElemBlocks; b++)
  {
    char elemType[MAX_LINE_LENGTH+1];
    int elsInBlock;
    int nodesPerEl;
    int numAttrs;
    int bid = blockIDs[b];

    ierr = ex_get_elem_block(exoID, bid, elemType, &elsInBlock,
      &nodesPerEl, &numAttrs);
    TEST_FOR_EXCEPT(ierr < 0);

    Array<int> connect(elsInBlock * nodesPerEl);

    ierr = ex_get_elem_conn(exoID, bid, &(connect[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    int n=0;

    for (int e=0; e<elsInBlock; e++, n+=nodesPerEl, count++)
    {
      if (dim==2)
      {
        mesh.addElement(elemGID[count], tuple(ptGID[connect[n]-1], ptGID[connect[n+1]-1], ptGID[connect[n+2]-1]), elemOwner[count], bid);
        SUNDANCE_VERB_HIGH("adding element=("
          << connect[n]-1 << ", " << connect[n+1]-1
          << ", " << connect[n+2]-1 << ")");
      }
      else
      {
        mesh.addElement(elemGID[count], 
          tuple(ptGID[connect[n]-1], ptGID[connect[n+1]-1], ptGID[connect[n+2]-1], ptGID[connect[n+3]-1]),
          elemOwner[count], bid);
      }
    }
  }


  /* Read the node sets */
  Array<int> nsIDs(numNodeSets);
  if (numNodeSets > 0)
  {
    ierr = ex_get_node_set_ids(exoID, &(nsIDs[0]));
  }
  TEST_FOR_EXCEPT(ierr < 0);
  for (int ns=0; ns<numNodeSets; ns++)
  {
    int nNodes;
    int nDist;
    int nsID = nsIDs[ns];
    ierr = ex_get_node_set_param(exoID, nsID, &nNodes, &nDist);
    TEST_FOR_EXCEPT(ierr < 0);
    Array<int> nodes(nNodes);
    ierr = ex_get_node_set(exoID, nsID, &(nodes[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    for (int n=0; n<nNodes; n++)
    {
      mesh.setLabel(0, nodes[n]-1, nsID);
    }
  }

    
  /* Read the side sets */
  Array<int> ssIDs(numSideSets);
  if (numSideSets > 0)
  {
    ierr = ex_get_side_set_ids(exoID, &(ssIDs[0]));
  }
  TEST_FOR_EXCEPT(ierr < 0);
  for (int ss=0; ss<numSideSets; ss++)
  {
    int nSides;
    int nDist;
    int ssID = ssIDs[ss];
    ierr = ex_get_side_set_param(exoID, ssID, &nSides, &nDist);
    TEST_FOR_EXCEPT(ierr < 0);
    Array<int> sides(nSides);
    Array<int> elems(nSides);
    ierr = ex_get_side_set(exoID, ssID, &(elems[0]), &(sides[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    for (int n=0; n<nSides; n++)
    {
      int elemID = elems[n];
      int facetNum = sides[n];
      int fsign;
      int sideLID = mesh.facetLID(dim, elemID-1, dim-1, facetNum-1,fsign);
      mesh.setLabel(dim-1, sideLID, ssID);
    }
  }


  /* Read the nodal attributes */
  int nNodalVars = 0;
  ierr = ex_get_var_param(exoID, "N", &nNodalVars);
  TEST_FOR_EXCEPT(ierr < 0);

  Array<Array<double> >& funcVals = *nodeAttributes();
  funcVals.resize(nNodalVars);

  for (int i=0; i<nNodalVars; i++)
  {
    int t = 1;
    funcVals[i].resize(mesh.numCells(0));
    ierr = ex_get_nodal_var(exoID, t, i+1, mesh.numCells(0), &(funcVals[i][0]));
    TEST_FOR_EXCEPT(ierr < 0);
  }

  /* Read the element attributes */
  int nElemVars = 0;
  ierr = ex_get_var_param(exoID, "E", &nElemVars);
  TEST_FOR_EXCEPT(ierr < 0);

  Array<Array<double> >& eFuncVals = *elemAttributes();
  eFuncVals.resize(nElemVars);

  for (int i=0; i<nElemVars; i++)
  {
    int t = 1;
    eFuncVals[i].resize(mesh.numCells(mesh.spatialDim()));
    ierr = ex_get_elem_var(exoID, t, i+1, 1, mesh.numCells(mesh.spatialDim()), &(eFuncVals[i][0]));
    TEST_FOR_EXCEPT(ierr < 0);
  }

  ierr = ex_close(exoID);
  TEST_FOR_EXCEPT(ierr < 0);

#endif
	return mesh;
}




void ExodusMeshReader::readParallelInfo(Array<int>& ptGID, 
  Array<int>& ptOwner,
  Array<int>& cellGID, 
  Array<int>& cellOwner) const
{
  int nPoints;
  int nElems;
  std::string line;
  Array<string> tokens;
  try
  {
    ptGID.resize(0);
    ptOwner.resize(0);
    cellGID.resize(0);
    cellOwner.resize(0);

    /* if we're running in parallel, read the info on processor 
     * distribution */
    if (nProc() > 1)
    {
      RCP<std::ifstream> parStream 
        = openFile(parFilename_, "parallel info");
     
      /* read the number of processors and the processor rank in 
       * the file. These must be consistent with the current number of
       * processors and the current rank */
      getNextLine(*parStream, line, tokens, '#');
      
      TEST_FOR_EXCEPTION(tokens.length() != 2, RuntimeError,
        "ExodusMeshReader::getMesh() expects 2 entries "
        "on the first line of .par file. In " 
        << parFilename_ << " I found \n[" << line << "]\n");

      int np = atoi(tokens[1]);
      int pid = atoi(tokens[0]);

      /* check consistency with the current number of
       * processors and the current rank */
      
      TEST_FOR_EXCEPTION(np != nProc(), RuntimeError,
        "ExodusMeshReader::getMesh() found "
        "a mismatch between the current number of processors="
        << nProc() << 
        "and the number of processors=" << np
        << "in the file " << parFilename_);

      TEST_FOR_EXCEPTION(pid != myRank(), RuntimeError,
        "ExodusMeshReader::getMesh() found "
        "a mismatch between the current processor rank="
        << myRank() << "and the processor rank="
        << pid << " in the file " << parFilename_);

      /* read the number of points */
      getNextLine(*parStream, line, tokens, '#');

      TEST_FOR_EXCEPTION(tokens.length() != 1, RuntimeError,
        "ExodusMeshReader::getMesh() requires 1 entry "
        "on the second line of .pxo file. Found line \n[" 
        << line << "]\n in file " << parFilename_);
      
      nPoints = StrUtils::atoi(tokens[0]);

      /* read the global ID and the owner PID for each point */
      ptGID.resize(nPoints);
      ptOwner.resize(nPoints);

      for (int i=0; i<nPoints; i++)
      {
        getNextLine(*parStream, line, tokens, '#');

        TEST_FOR_EXCEPTION(tokens.length() != 3, RuntimeError,
          "ExodusMeshReader::getMesh() requires 3 "
          "entries on each line of the point section in "
          "the .pxo file. Found line \n[" << line
          << "]\n in file " << parFilename_);

        ptGID[i] = StrUtils::atoi(tokens[1]);
        ptOwner[i] = StrUtils::atoi(tokens[2]);
      }


      /* Read the number of elements */

      getNextLine(*parStream, line, tokens, '#');

      TEST_FOR_EXCEPTION(tokens.length() != 1, RuntimeError,
        "ExodusMeshReader::getMesh() requires 1 entry "
        "on the cell count line of .pxo file. Found line \n[" 
        << line << "]\n in file " << parFilename_);

      nElems = StrUtils::atoi(tokens[0]);

      SUNDANCE_OUT(this->verb() > 1,
        "read nElems = " << nElems);


      /* read the global ID and the owner PID for each element */

      cellGID.resize(nElems);
      cellOwner.resize(nElems);
      for (int i=0; i<nElems; i++)
      {
        getNextLine(*parStream, line, tokens, '#');

        TEST_FOR_EXCEPTION(tokens.length() != 3, RuntimeError,
          "ExodusMeshReader::getMesh() requires 3 "
          "entries on each line of the element section in "
          "the .pxo file. Found line \n[" << line
          << "]\n in file " << parFilename_);

        cellGID[i] = StrUtils::atoi(tokens[1]);
        cellOwner[i] = StrUtils::atoi(tokens[2]);
      }
    }

    nPoints = ptGID.length();
    nElems = cellGID.length();
  }
  catch(std::exception& e)
  {
    SUNDANCE_TRACE(e);
  }
}
