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

#ifndef SUNDANCE_TRIANGLEMESHREADER_H
#define SUNDANCE_TRIANGLEMESHREADER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"

namespace Sundance
{
  using namespace Teuchos;
  
  /**
   * TriangleMeshReader reads a mesh stored in Shewchuk's Triangle format.
   * This format is documented at 
   * <A HREF="http://www-2.cs.cmu.edu/~quake/triangle.html"> 
   * the Triangle homepage. 
   * </A>
   * This reader expects to find node information in <tt>.node</tt> files
   * and element information in <tt>.ele</tt> files. The <tt> filename </tt>
   * constructor argument is the stem of the filenames, and so that 
   * a reader constructed with filename <tt>joe</tt> will look for node and
   * element data in <tt>joe.node</tt> and <tt>joe.ele</tt> respectively.
   * Node and element
   * attributes are read if present, and can be accessed with the 
   * <tt>getAttributes()</tt> method of <tt>MeshSource.</tt>
   * 
   * <h4> Parallel extensions </h4>
   * We have extended the Triangle format to deal with distributed meshes.
   * A TriangleMeshReader is constructed with an MPIComm object, and if
   * that communicator has more than one processor the mesh is assumed
   * to be split into files, one for each processor. Data
   * on mesh "joe" for the <i>nnn</i>-th processor will be found in the files
   * <ul>
   * <li> <tt>joe.node.</tt><i>nnn</i>
   * <li> <tt>joe.ele.</tt><i>nnn</i>
   * <li> <tt>joe.par.</tt><i>nnn</i>
   * </ul>
   * The <tt>.node.</tt><i>nnn</i> and <tt>.ele.</tt><i>nnn</i> files contain the
   * node and element data for the part of the mesh seen 
   * by the <i>nnn</i>-th
   * processor. The node and element 
   * numberings given in those two files are <b>local</b> indexes.
   * The <tt>.par.</tt><i>nnn</i> contains node and element 
   * ownership information for the part of the mesh seen 
   * by the <i>nnn</i>-th
   * processor. 
   *
   * <br> 
   *
   * A <tt>.par</tt> file is formatted as follows:
   * <ul>
   * <li> First line: <tt> rank numProcs </tt>
   * <li> Second line: <tt> numPoints </tt>
   * <li> Next <i> nPoints </i> lines: <tt> ptLID ptGID ptOwnerRank </tt>
   * <li> Next line: <tt> numElems </tt>
   * <li> Next <i> nElems </i> lines: <tt> elemLID elemGID elemOwnerRank </tt>
   * </ul>
   * 
   */
  class TriangleMeshReader : public MeshReaderBase
  {
  public:
    /** */
    TriangleMeshReader(const std::string& filename, 
                       const MeshType& meshType,
                       const MPIComm& comm = MPIComm::world());

    /** Construct from a ParameterList */
    TriangleMeshReader(const ParameterList& params);

    /** virtual dtor */
    virtual ~TriangleMeshReader(){;}


    /** Create a mesh */
    virtual Mesh fillMesh() const ;

    /** Print a short descriptive std::string */
    virtual std::string description() const 
    {return "TriangleMeshReader[file=" + filename() + "]";}
      

    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

  private:
    /** */
    void readParallelInfo(Array<int>& ptGID, Array<int>& ptOwner,
                          Array<int>& elemGID, Array<int>& elemOwner) const ;

    /** */
    Mesh readNodes(Array<int>& ptGID,
                   Array<int>& ptOwner) const ;

    void readSides(Mesh& mesh) const ;

    /** */
    void readElems(Mesh& mesh,
                   const Array<int>& nodeGID,
                   Array<int>& elemGID,
                   Array<int>& elemOwner) const ;
    

    /** */
    std::string nodeFilename_;

    /** */
    std::string elemFilename_;

    /** */
    std::string parFilename_;

    /** */
    std::string sideFilename_;

    /** */
    mutable int offset_;
  };
}

#endif
