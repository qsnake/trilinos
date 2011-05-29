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

#ifndef SUNDANCE_FILEIO_CHACOPARTITIONER_H
#define SUNDANCE_FILEIO_CHACOPARTITIONER_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceSerialPartitionerBase.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{
  /**
   * This partitioner makes a system call to run Chaco, passing partitioning
   * information to it through file IO.
   */
class FileIOChacoPartitioner : public SerialPartitionerBase
  {
  public:

    /** Construct an empty mesh filter object */
    FileIOChacoPartitioner(const std::string& filename);



    /** */
    virtual void getAssignments(const Mesh& mesh, int np, 
      Array<int>& assignments) const ;

    /** */
    void writeGraph(const Mesh& mesh) const ;
    
    /** */
    void runChaco(int np) const ;

  private:

  /** Determine whether a line is empty */
  bool isEmptyLine(const std::string& x) const ;

  /** 
   * Read the next non-empty, non-comment line from a stream
   * @param is the stream from which to get the line
   * @param line upon return, filled in with the line that was read
   * @param tokens array of space-separated tokens in the line
   * @param comment a character indicating that everything after it
   * is a comment
   */
  bool getNextLine(std::istream& is, std::string& line,
    Array<string>& tokens,
    char comment) const ;

    std::string filename_;
  };
}

#endif
