#include "SundanceFileIOChacoPartitioner.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

using std::ofstream;
using std::ifstream;
using std::endl;


FileIOChacoPartitioner::FileIOChacoPartitioner(const std::string& filename)
  : filename_(filename)
{}

void FileIOChacoPartitioner::writeGraph(const Mesh& mesh) const 
{
  Array<Array<int> > neighbors;
  int nEdges;

  getNeighbors(mesh, neighbors, nEdges);

  std::string gf = filename_ + ".graph";
  ofstream os(gf.c_str());

  os << neighbors.size() << " " << nEdges << std::endl;

  for (int i=0; i<neighbors.size(); i++)
  {
    for (int j=0; j<neighbors[i].size(); j++) 
    {
      if (j > 0) os << " ";
      os << neighbors[i][j]+1; // need unit offset here for Chaco
    }
    os << "\n";
  }
}


void FileIOChacoPartitioner::runChaco(int np) const 
{
  ofstream pf("User_Params");
  pf << 
    "OUTPUT_ASSIGN=true\n"
    "PROMPT=false\n"
    "ARCHITECTURE=1\n"
    "REFINE_PARTITION=4\n"
    "REFINE_MAP=true\n"
    "KL_BAD_MOVES=20\n"
    "KL_NTRIES_BAD=10\n"
    "KL_IMBALANCE=0.02\n"
    "INTERNAL_VERTICES=true\n"
    "MATCH_TYPE=2\n"
    "HEAVY_MATCH=true\n"
    "TERM_PROP=true\n"
    "COARSE_NLEVEL_KL=1\n"
    "COARSEN_RATIO_MIN=0.7\n"
    "CUT_TO_HOP_COST=1.0\n"
    "RANDOM_SEED=12345\n" << std::endl;

  ofstream chIn("chacoInput");
  chIn << filename_ + ".graph\n" << filename_ + ".assign\n1\n100\n"
       << np << "\n1\nn" << std::endl;

  int status = system("chaco < chacoInput");
  TEST_FOR_EXCEPTION(status < 0, RuntimeError, "error detected in system call to run chaco");
}

bool FileIOChacoPartitioner::isEmptyLine(const std::string& x) const 
{
  return x.length()==0 || StrUtils::isWhite(x);
}

bool FileIOChacoPartitioner::getNextLine(std::istream& is, std::string& line,
                                         Array<string>& tokens,
                                         char comment) const 
{
  bool rtn = false;
  while ((rtn=StrUtils::readLine(is, line)))
    {
      if (line.length() > 0) line = StrUtils::before(line,comment);
      if (isEmptyLine(line)) continue;
      if (line.length() > 0) break;
    }
  tokens = StrUtils::stringTokenizer(line);
  return rtn;
}

void FileIOChacoPartitioner::getAssignments(const Mesh& mesh, int np,
  Array<int>& assignments) const 
{
  writeGraph(mesh);
  runChaco(np);

  std::string af = filename_ + ".assign";
  ifstream is(af.c_str());

  std::string line;
  Array<string> tokens;
    
  while (getNextLine(is, line, tokens, '#'))
  {
    assignments.append(StrUtils::atoi(tokens[0]));
  }
}

