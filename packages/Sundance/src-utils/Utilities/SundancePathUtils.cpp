/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundancePathUtils.hpp"
#include "SundanceDefaultPath.hpp"
#include <unistd.h>
#ifndef _MSC_VER
#include <sys/unistd.h>
#endif

using Teuchos::Array;

using std::ifstream;

namespace Sundance
{
string searchForFile(const std::string& name)
{
  std::string pathSep = "/";
  Array<string> path = parsePathStr();

  if (name.length() && name[0]=='/') return name; // Use absolute path!
  for (int i=0; i<path.size(); i++)
  {
    ifstream fileToTry((path[i] + pathSep + name).c_str());
    if (!fileToTry) continue;
    return path[i] + pathSep + name;
  }

  TEST_FOR_EXCEPTION(true, std::runtime_error, "could not find file "
    << name << " in path " << path);
}

string getPathStr() 
{
  char* pathEnvStr = getenv("SUNDANCE_PATH");
  char* pyPathEnvStr = getenv("PYTHONPATH");
  std::string path;
  
  if (pathEnvStr == NULL) 
  {
    path = defaultSundancePath();
  }
  else
  {
    path = pathEnvStr;
  }
  if (pyPathEnvStr!=NULL)
  {
    path = std::string(pyPathEnvStr) + ":" + path; 
  }
  return path;
}

Array<string> parsePathStr() 
{
  std::string pathStr = getPathStr();
  
  Array<string> rtn;

  unsigned int begin;
  unsigned int end;
  
  begin = pathStr.find_first_not_of(":");
  
  while (begin < pathStr.length())
  {
    end = pathStr.find_first_of(":", begin);

    rtn.append(pathStr.substr(begin, end-begin));
    begin = pathStr.find_first_not_of(":", end);
  }

  return rtn;
}
}
