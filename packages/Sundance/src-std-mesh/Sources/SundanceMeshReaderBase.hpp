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

#ifndef SUNDANCE_MESHREADERBASE_H
#define SUNDANCE_MESHREADERBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "Teuchos_StrUtils.hpp"

namespace Sundance
{

/**
 * MeshReaderBase is a base class for mesh sources that get a mesh
 * from a file. It provides several utilities for parsing lines
 * from mesh files. 
 */
class MeshReaderBase : public MeshSourceBase
{
public:
  /** Construct with a filename */
  MeshReaderBase(const std::string& filename,
    const MeshType& meshType,
    const MPIComm& comm)
    : MeshSourceBase(meshType, comm), filename_(filename)
    {}

  /** Construct from a parameter list */
  MeshReaderBase(const ParameterList& params);

  /** */
  virtual ~MeshReaderBase(){;}

protected:
  /** access to the filename */
  const std::string& filename() const {return filename_;}

  /** convert a std::string to its integer value */
  int atoi(const std::string& x) const ;

  /** convert a std::string to its double value */
  double atof(const std::string& x) const ;

  /** Determine whether a line is empty */
  bool isEmptyLine(const std::string& x) const ;

  /** Open a file "fname" and check for success.
   * @param fname name of the file to be opened
   * @param description a description of the file, e.g., "node file",
   * to be included in any error messages generated.  
   **/
  RCP<std::ifstream> openFile(const std::string& fname, 
    const std::string& description) const ;

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
private:
  std::string filename_;
  mutable int nVertexVars_;
  mutable Array<double> vertexVars_;
};
}



#endif
