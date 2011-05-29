#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


ExodusNetCDFMeshReader::ExodusNetCDFMeshReader(const std::string& fname,
                                               const MeshType& meshType,
                                               const MPIComm& comm)
  : MeshReaderBase(fname, meshType, comm)
{
  TEST_FOR_EXCEPTION(nProc() > 1, RuntimeError,
                     "ExodusNetCDFMeshReader not useable with parallel meshes");
}

ExodusNetCDFMeshReader::ExodusNetCDFMeshReader(const ParameterList& params)
  : MeshReaderBase(params)
{
  TEST_FOR_EXCEPTION(nProc() > 1, RuntimeError,
                     "ExodusNetCDFMeshReader not useable with parallel meshes");
}




Mesh ExodusNetCDFMeshReader::fillMesh() const 
{
  Mesh mesh;

  RCP<std::ifstream> is = openFile(filename(), "NetCDF");

  std::string line;
  Array<string> tokens;

  /* read the header line */
  getNextLine(*is, line, tokens, '#');

  TEST_FOR_EXCEPTION(tokens[0] != "netcdf", RuntimeError,
                     "ExodusNetCDF reader expected to find [netcdf] as first "
                     "token, found " << tokens[0]);

  /* read the list of dimensions */
  getNextLine(*is, line, tokens, '#');


  TEST_FOR_EXCEPTION(tokens[0] != "dimensions:", RuntimeError,
                     "ExodusNetCDF reader expected to find [dimension:] as first "
                     "token, found " << tokens[0]);
  
  int nElem = 0;
  int nNodes = 0;
  int nElemBlocks = 0;
  int nSideSets = 0;
  int nNodeSets = 0;
  int dimension = 0 ;
  Array<int> blockSizes;
  Array<int> sideSetSizes;
  Array<int> nodeSetSizes;
  Array<int> nodesPerElem;
  while (true)
    {
      getNextLine(*is, line, tokens, '#');

      if (tokens[0] == "variables:") break;
      std::string keyword = tokens[0];
      std::string equals = tokens[1];
      TEST_FOR_EXCEPTION(equals!="=", RuntimeError, "ExodusNetCDF reader "
                         "expected [=] as second token, found "
                         << equals);

      TEST_FOR_EXCEPTION(tokens.size() < 4, RuntimeError,
                         "ExodusNetCDF reader found a dimension line with "
                         "fewer than 4 tokens");

      int val = atoi(tokens[2]);
      if (keyword=="num_dim")
        {
          dimension = val;
        }
      else if (keyword=="num_nodes")
        {
          nNodes = val;
        }
      else if (keyword=="num_elem")
        {
          nElem = val;
        } 
      else if (keyword=="num_el_blk")
        {
          nElemBlocks = val;
          blockSizes.resize(nElemBlocks);
          nodesPerElem.resize(nElemBlocks);
        }
      else if (keyword=="num_side_sets")
        {
          nSideSets = val;
          sideSetSizes.resize(nSideSets);
          std::cerr << "num side sets = " << nSideSets << std::endl;
        }
      else if (keyword=="num_node_sets")
        {
          nNodeSets = val;
          nodeSetSizes.resize(nNodeSets);
          std::cerr << "num node sets = " << nNodeSets << std::endl;
        }
      else
        {
          for (int i=0; i<nElemBlocks; i++)
            {
              if (keyword=="num_el_in_blk" + Teuchos::toString(i+1))
                {
                  blockSizes[i] = val;
                  break;
                }
              if (keyword=="num_nod_per_el" + Teuchos::toString(i+1))
                {
                  nodesPerElem[i] = val;
                  break;
                }
            }
          for (int i=0; i<nSideSets; i++)
            {
              if (keyword=="num_side_ss" + Teuchos::toString(i+1))
                {
                  sideSetSizes[i] = val;
                  break;
                }
            }
          for (int i=0; i<nNodeSets; i++)
            {
              if (keyword=="num_nod_ns" + Teuchos::toString(i+1))
                {
                  nodeSetSizes[i] = val;
                  break;
                }
            }
        }
    }

  

  /* skip until we find [data:] */
  while (true)
    {
      getNextLine(*is, line, tokens, '#');

      if (tokens[0] == "data:") break;
    }

  /* read the data */
  

  Array<double> coords;
  coords.reserve(nNodes*dimension);
  Array<Array<int> > connect(nElemBlocks);
  Array<Array<int> > sideSetElems(nSideSets);
  Array<Array<int> > sideSetFacets(nSideSets);
  Array<Array<int> > nodeSetNodes(nNodeSets);

  bool doneWithData = false;
  while(!doneWithData)
    {
      if (!getNextLine(*is, line, tokens, '#')) 
        {
          doneWithData = true;
          break;
        }

      if (tokens.size()==0) 
        {

          doneWithData=true;
          break;
        }
      
      if (tokens[0]=="coord")
        {
          bool done = false;
          for (int i=1; i<tokens.size(); i++)
            {
              if (tokens[i] == "=") continue;
              if (tokens[i] == ";") 
                {
                  done = true;
                  break;
                }
              coords.append(atof(StrUtils::before(tokens[i], ",")));
            }
          while (!done)
            {
              if (!getNextLine(*is, line, tokens, '#'))
                {
                  doneWithData = true;
                  break;
                }

              for (int i=0; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  coords.append(atof(StrUtils::before(tokens[i], ",")));
                }
            }
          continue;
        }
      /* see if the line lists connectivity data for a block */
      for (int b=0; b<nElemBlocks; b++)
        {
          if (tokens[0] == "connect" + Teuchos::toString(b+1))
            {
              connect[b].reserve(blockSizes[b]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  connect[b].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }

                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      connect[b].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              continue;
            }
        }

      /* see if the line lists side set element numbers */
      for (int s=0; s<nSideSets; s++)
        {
          if (tokens[0] == "elem_ss" + Teuchos::toString(s+1))
            {
              sideSetElems[s].reserve(sideSetSizes[s]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  sideSetElems[s].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }

                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      sideSetElems[s].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              continue;
            }
        }

      /* see if the line lists side set facet numbers */
      for (int s=0; s<nSideSets; s++)
        {
          if (tokens[0] == "side_ss" + Teuchos::toString(s+1))
            {
              sideSetFacets[s].reserve(sideSetSizes[s]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  sideSetFacets[s].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }
                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      sideSetFacets[s].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              continue;
            }
        }

      /* see if the line lists node set node numbers */
      for (int s=0; s<nNodeSets; s++)
        {
          if (tokens[0] == "node_ns" + Teuchos::toString(s+1))
            {
              nodeSetNodes[s].reserve(nodeSetSizes[s]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  nodeSetNodes[s].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }
                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      nodeSetNodes[s].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              continue;
            }
        }
    }



  /* now we can build the mesh */
  mesh = createMesh(dimension);

  /* do some consistency checks */
  TEST_FOR_EXCEPTION(dimension * nNodes != coords.size(), RuntimeError,
                     "bad coordinate array in exodus reader");

  for (int b=0; b<nElemBlocks; b++)
    {
      TEST_FOR_EXCEPTION(blockSizes[b]*nodesPerElem[b] != connect[b].size(), RuntimeError,
                         "bad connectivity array for block " << b << " in exodus reader");
    }

  for (int s=0; s<nSideSets; s++)
    {
      TEST_FOR_EXCEPTION(sideSetElems[s].size() != sideSetFacets[s].size(), RuntimeError,
                         "inconsistent side set specification for ss=" << s 
                         << " in exodus reader ");
    }

  /* add the points to the mesh */
  for (int n=0; n<nNodes; n++)
    {
      Point x;
      if (dimension==2)
        {
          x = Point(coords[n], coords[nNodes+n]);
        }
      else
        {
          x = Point(coords[n], coords[nNodes+n], coords[2*nNodes+n]);
        }
      mesh.addVertex(n, x, 0, 0);
    }


  /* add the elements */
  int lid=0;
  for (int b=0; b<nElemBlocks; b++)
    {
      int n = 0;
      for (int e=0; e<blockSizes[b]; e++, n+=nodesPerElem[b], lid++)
        {
          if (dimension==2)
            {
              mesh.addElement(lid, tuple(connect[b][n]-1, connect[b][n+1]-1, connect[b][n+2]-1), 0, b+1);
            }
          else
            {
              mesh.addElement(lid, 
                              tuple(connect[b][n]-1, connect[b][n+1]-1, connect[b][n+2]-1, connect[b][n+3]-1),
                              0, b+1);
            }
        }
    }
  
  /* label the side sets */

  for (int s=0; s<nSideSets; s++)
    {
      for (int n=0; n<sideSetSizes[s]; n++)
        {
          int elemID = sideSetElems[s][n];
          int facetNum = sideSetFacets[s][n];
          int fsign;
          int sideLID = mesh.facetLID(dimension, elemID-1, dimension-1, facetNum-1,fsign);
          mesh.setLabel(dimension-1, sideLID, s+1);
        }
    }

  /* label the node sets */
  for (int s=0; s<nNodeSets; s++)
    {
      for (int n=0; n<nodeSetSizes[s]; n++)
        {
          int nodeNum = nodeSetNodes[s][n]-1;
          mesh.setLabel(0, nodeNum, s+1);
        }
    }

	return mesh;
}

