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

#ifndef SUNDANCE_OUT_H
#define SUNDANCE_OUT_H

#include "SundanceDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_MPIComm.hpp"

namespace Teuchos
{
template <class T> class Array;
}
namespace Sundance
{
class Tabs;
using namespace Teuchos;

/**
 *
 */
class Out
{
public:
      
  static void println(const std::string& str) 
    {
      if (hasLogFile()) *logFile() << str << std::endl;
      if (!suppressStdout()) os() << str << std::endl;
    }

  static void setLogFile(const std::string& filename)
    {
      logFile() = rcp(new std::ofstream(filename.c_str()));
      hasLogFile() = true;
    }

  static FancyOStream& os()
    {
      static RCP<std::ostream> os = rcp(&std::cout, false);
      static RCP<FancyOStream> rtn = fancyOStream(os);
      static bool first = true;
      if (first)
      {
        rtn->setShowProcRank(true);
        first = false;
      }
      return *rtn;
    }

  static FancyOStream& root()
    {
      static bool isRoot = MPIComm::world().getRank()==0;
      static RCP<FancyOStream> blackHole
        = rcp(new FancyOStream(rcp(new oblackholestream())));

      if (isRoot)
      {
        return os();
      }
      else
      {
        return *blackHole;
      }
    }

  static bool& suppressStdout() {static bool rtn=false; return rtn;}

private:
  static bool& hasLogFile() {static bool rtn=false; return rtn;}
  static RCP<std::ostream>& logFile() {static RCP<std::ostream> rtn; return rtn;}
      
};


/** */
void writeTable(std::ostream& os, const Tabs& tab,
  const Array<double>& a, int cols);
/** */
void writeTable(std::ostream& os, const Tabs& tab, 
  const Array<int>& a, int cols);

}

#define SUNDANCE_OUT(test, msg) \
  { \
    if (test) \
    { \
      TeuchosOStringStream omsg; \
      omsg << msg; \
      Out::println(omsg.str());                 \
    } \
  }


#define SUNDANCE_VERB_EXTREME(msg) SUNDANCE_MSG4(this->verb(), msg)
#define SUNDANCE_VERB_HIGH(msg) SUNDANCE_MSG3(this->verb(), msg)
#define SUNDANCE_VERB_MEDIUM(msg) SUNDANCE_MSG2(this->verb(), msg)
#define SUNDANCE_VERB_LOW(msg) SUNDANCE_MSG1(this->verb(), msg)

#define SUNDANCE_HEADER_LINE "\n------------------------------------------------------------------\n"

#define SUNDANCE_MSG(context, level, msg) SUNDANCE_OUT(this->verbLevel(context) >= level, msg)

#define SUNDANCE_LEVEL1(context, msg) SUNDANCE_MSG(context, 1, msg)

#define SUNDANCE_LEVEL2(context, msg) SUNDANCE_MSG(context, 2, msg)

#define SUNDANCE_LEVEL3(context, msg) SUNDANCE_MSG(context, 3, msg)

#define SUNDANCE_LEVEL4(context, msg) SUNDANCE_MSG(context, 4, msg)

#define SUNDANCE_LEVEL5(context, msg) SUNDANCE_MSG(context, 5, msg)


#define SUNDANCE_MSG1(level, msg) SUNDANCE_OUT(level >= 1, msg)

#define SUNDANCE_MSG2(level, msg) SUNDANCE_OUT(level >= 2, msg)

#define SUNDANCE_MSG3(level, msg) SUNDANCE_OUT(level >= 3, msg)

#define SUNDANCE_MSG4(level, msg) SUNDANCE_OUT(level >= 4, msg)

#define SUNDANCE_MSG5(level, msg) SUNDANCE_OUT(level >= 5, msg)

#define SUNDANCE_BANNER1(level, tab, msg) \
  SUNDANCE_MSG1(level, tab \
    << "===================================================================");\
  SUNDANCE_MSG1(level, tab << std::endl << tab \
    << "  " << msg); \
  SUNDANCE_MSG1(level, tab << std::endl << tab\
    << "===================================================================");


#define SUNDANCE_BANNER2(level, tab, msg) \
  SUNDANCE_MSG2(level, tab \
    << "-------------------------------------------------------------------");\
  SUNDANCE_MSG2(level, tab << msg); \
  SUNDANCE_MSG2(level, tab\
    << "-------------------------------------------------------------------");



#define SUNDANCE_BANNER3(level, tab, msg) \
  SUNDANCE_MSG3(level, tab \
    << "-------------------------------------------------------------------");\
  SUNDANCE_MSG3(level, tab << std::endl << tab \
    << msg); \
  SUNDANCE_MSG3(level, tab << std::endl << tab\
    << "-------------------------------------------------------------------");

#endif
