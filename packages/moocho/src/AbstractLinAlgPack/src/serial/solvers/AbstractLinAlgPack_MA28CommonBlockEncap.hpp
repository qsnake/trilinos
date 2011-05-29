// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
//
// These classes are used to encapsulate the copy of ma28 common block data from
// on set to another.

#ifndef MA28_COMMON_BLOCK_ENCAP_H
#define MA28_COMMON_BLOCK_ENCAP_H

#include <iostream>

#include "AbstractLinAlgPack_MA28_CppDecl.hpp"

namespace MA28_Cpp {

  using MA28_CppDecl::MA28ED_struct;
  using MA28_CppDecl::MA28FD_struct;
  using MA28_CppDecl::MA28GD_struct;
  using MA28_CppDecl::MA28HD_struct;
  using MA28_CppDecl::MA30ED_struct;
  using MA28_CppDecl::MA30FD_struct;
  using MA28_CppDecl::MA30GD_struct;
  using MA28_CppDecl::MA30HD_struct;
  using MA28_CppDecl::MA30ID_struct;
  using MA28_CppDecl::MC23BD_struct;

class MA28CommonBlockStorage; // forward declaration

// This is a class for storing the references to MA28 common blocks.
// It is used to simplify copy of common block members.
class MA28CommonBlockReferences {
public:
  // Construct with references.
  MA28CommonBlockReferences(
    MA28ED_struct& ma28ed,
    MA28FD_struct& ma28fd,
    MA28GD_struct& ma28gd,
    MA28HD_struct& ma28hd,
    MA30ED_struct& ma30ed,
    MA30FD_struct& ma30fd,
    MA30GD_struct& ma30gd,
    MA30HD_struct& ma30hd,
    MA30ID_struct& ma30id,
    MC23BD_struct& mc23bd
    ) :
    ma28ed_(ma28ed),
    ma28fd_(ma28fd),
    ma28gd_(ma28gd),
    ma28hd_(ma28hd),
    ma30ed_(ma30ed),
    ma30fd_(ma30fd),
    ma30gd_(ma30gd),
    ma30hd_(ma30hd),
    ma30id_(ma30id),
    mc23bd_(mc23bd)
  {}

  // Use default copy constructor ....

  // public reference members
  MA28ED_struct& ma28ed_;
  MA28FD_struct& ma28fd_;
  MA28GD_struct& ma28gd_;
  MA28HD_struct& ma28hd_;
  MA30ED_struct& ma30ed_;
  MA30FD_struct& ma30fd_;
  MA30GD_struct& ma30gd_;
  MA30HD_struct& ma30hd_;
  MA30ID_struct& ma30id_;
  MC23BD_struct& mc23bd_;

  // Assignment operator
  MA28CommonBlockReferences& operator=(const MA28CommonBlockStorage& ma28cbs);

  // Use default assignment operator ....

  // Output common block info to a ostream
  void dump_values(std::ostream& o);
private:
  MA28CommonBlockReferences(); // not defined and not to be called.
};

// This is a storage class for MA28 common blocks.  It is ment to be
// used by objects that wish to save the state of MA28.  The default
// C++ copy constructor and assignment operator are use since a 
// memberwise copy of member data is exactly what I want.
class MA28CommonBlockStorage {
public:
  // Use default constructor ....
  // Use default copy constructor .....
  
  // Copy the values
  MA28CommonBlockStorage(const MA28CommonBlockReferences& cb) :
    ma28ed_(cb.ma28ed_),
    ma28fd_(cb.ma28fd_),
    ma28gd_(cb.ma28gd_),
    ma28hd_(cb.ma28hd_),
    ma30ed_(cb.ma30ed_),
    ma30fd_(cb.ma30fd_),
    ma30gd_(cb.ma30gd_),
    ma30hd_(cb.ma30hd_),
    ma30id_(cb.ma30id_),
    mc23bd_(cb.mc23bd_)
  {}

  // Public storage members
  MA28ED_struct ma28ed_;
  MA28FD_struct ma28fd_;
  MA28GD_struct ma28gd_;
  MA28HD_struct ma28hd_;
  MA30ED_struct ma30ed_;
  MA30FD_struct ma30fd_;
  MA30GD_struct ma30gd_;
  MA30HD_struct ma30hd_;
  MA30ID_struct ma30id_;
  MC23BD_struct mc23bd_;

  // Assignment operator
  MA28CommonBlockStorage& operator=(const MA28CommonBlockReferences& ma28cbs);

  // Use default assignment operator

  // Output common block info to a ostream
  void dump_values(std::ostream& o);

};

}	// end namespace MA28_Cpp

#endif // MA28_COMMON_BLOCK_ENCAP_H
