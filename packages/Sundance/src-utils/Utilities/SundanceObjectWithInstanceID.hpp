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

#ifndef SUNDANCE_OBJECTWITHINSTANCEID_H
#define SUNDANCE_OBJECTWITHINSTANCEID_H

#include "SundanceDefs.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
    /**
     * ObjectWithInstanceID provides a common method for the
     * generation of instance-specific ID numbers. Subclasses will inherit
     * the id() method, and instances of those subclasses will be
     * given a unique ID at construction time.
     * 
     * <h4> Design note: </h4> By templating on the derived type, 
     * we can give each derived type its own sequence of ID numbers.
     */
    template <class T>
    class ObjectWithInstanceID
    {
    public:
      /** Empty ctor will assign ID at construction time */
      ObjectWithInstanceID() : id_(nextID()) {;} 

      /** Return this object's ID number */
      int id() const {return id_;}

    private:
      /** Generate the next ID in the sequence */
      static int& nextID() {static int rtn=0; rtn++; return rtn;}

      /** */
      int id_;
    };

}


#endif  /* DOXYGEN_DEVELOPER_ONLY */   
#endif
