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

#ifndef SUNDANCE_FIELDWRITER_H
#define SUNDANCE_FIELDWRITER_H

#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{
  /**
   * FieldWriter is the user level object for writing fields and meshes
   * to output stream. 
   *
   * <h4> Example: </h4> Write fields u0 and w0 to a VTK file "results.vtu"
   * \code
   * FieldWriter vtkWriter = new VTKWriter("results");
   * vtkWriter.addField(u0);
   * vtkWriter.addField(w0);
   * vtkWriter.write();
   * \endcode
   *
   * <h4> Example: </h4> Write verbose mesh information to cout
   * \code
   * FieldWriter writer = new VerboseFieldWriter();
   * writer.addMesh(mesh);
   * writer.write();
   * \endcode
   */
  class FieldWriter : public Sundance::Handle<FieldWriterBase>
  {
  public:
    /* Boilerplate handle ctors */
    HANDLE_CTORS(FieldWriter, FieldWriterBase);

    /** add a mesh to the list of things to be written */
    void addMesh(const Mesh& mesh) const ;

    /** add a field, tagging it with the given std::string as a name */
    void addField(const std::string& name, 
                  const Handle<FieldBase>& field) ;

    /** set the numerical value to be written at cells on which
     * a field is undefined. Default value is 0.0. */
    void setUndefinedValue(const double& x);

    /** 
     * Tell the writer to pretend that it is running as one of nProc processes
     * with the given rank. This is used when partitioning meshes in serial.
     */
    void impersonateParallelProc(int nProc, int rank)
      {
        ptr()->impersonateParallelProc(nProc, rank);
      }
    
    

    /** write to stream */
    void write() const ;
  };
}

#endif
