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

#ifndef SUNDANCE_EXODUSMESHREADER_H
#define SUNDANCE_EXODUSMESHREADER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
  using namespace Teuchos;
  
  /**
   * ExodusMeshReader reads a mesh from an ExodusII file.
   */
  class ExodusMeshReader : public MeshReaderBase
  {
  public:
    /** */
    ExodusMeshReader(const std::string& filename, 
                           const MeshType& meshType,
                           const MPIComm& comm = MPIComm::world());

    /** virtual dtor */
    virtual ~ExodusMeshReader(){;}


    /** Create a mesh */
    virtual Mesh fillMesh() const ;

    /** Print a short descriptive std::string */
    virtual std::string description() const 
    {return "ExodusMeshReader[file=" + filename() + "]";}

    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

  private:
    /** */
    void readParallelInfo(Array<int>& ptGID, Array<int>& ptOwner,
                          Array<int>& elemGID, Array<int>& elemOwner) const ;
    
    /** */
    std::string exoFilename_;
    /** */
    std::string parFilename_;
    
                        };
}

#endif
