// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_HYBRIDPLATFORM_HPP
#define TPETRA_HYBRIDPLATFORM_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterList.hpp>
#include <string>
#include <cstdio> // for std::sscanf

#include <Kokkos_SerialNode.hpp>
#ifdef HAVE_KOKKOS_TBB
#include <Kokkos_TBBNode.hpp>
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include <Kokkos_TPINode.hpp>
#endif
#ifdef HAVE_KOKKOS_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
#endif

namespace Tpetra {

	//! \brief A platform class for hybrid nodes.
  /*!
    This class is templated on two types, those of the two underlying Nodes.
    In this way, the HybridPlatform is compiled with support for a particular 
    hybrid architecture.
   */
  class HybridPlatform : public Teuchos::Describable {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      HybridPlatform(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, Teuchos::ParameterList &pl);

      //! Destructor
      ~HybridPlatform();

      //@}

      //! @name Class Query, Creation and Accessor Methods
      //@{ 

      //! Comm Instance
      Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      void createNode();

      //! Run user code with the runtime-selected Node type.
      template <template <class Node> class UserCode> 
      void runUserCode();

      //@}

    private:
      HybridPlatform(const HybridPlatform &platform); // not supported
      const Teuchos::RCP<const Teuchos::Comm<int> > comm_;
      Teuchos::ParameterList instList_;
      Teuchos::RCP<Kokkos::SerialNode>    serialNode_;
      bool nodeCreated_;
#ifdef HAVE_KOKKOS_TBB
      Teuchos::RCP<Kokkos::TBBNode>       tbbNode_;
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
      Teuchos::RCP<Kokkos::TPINode>       tpiNode_;
#endif
#ifdef HAVE_KOKKOS_THRUST
      Teuchos::RCP<Kokkos::ThrustGPUNode> thrustNode_;
#endif

      enum NodeType {
        SERIALNODE
#ifdef HAVE_KOKKOS_TBB
        , TBBNODE
#endif        
#ifdef HAVE_KOKKOS_THREADPOOL
        , TPINODE
#endif        
#ifdef HAVE_KOKKOS_THRUST
        , THRUSTGPUNODE
#endif        
      } nodeType_;
  };

  HybridPlatform::HybridPlatform(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, Teuchos::ParameterList &pl)
  : comm_(comm)
  , nodeCreated_(false)
  , nodeType_(SERIALNODE)
  {
    // ParameterList format:
    // 
    // Node designation sublists have a name beginning with one of the following: % = [
    // and satisfying the following format:
    //   %M=N    is satisfied if mod(myrank,M) == N
    //   =N      is satisfied if myrank == N
    //   [M,N]   is satisfied if myrank \in [M,N]
    // 
    // A node designation sublist must have a parameter entry of type std::string named "NodeType". The value indicates the type of the Node.
    // The activated node designation sublist will be passed to the Node constructor.
    // 
    // For example:
    // "%2=0"  ->  
    //    NodeType     = "Kokkos::ThrustGPUNode"
    //    DeviceNumber = 0
    //    Verbose      = 1
    // "%2=1"  ->
    //    NodeType     = "Kokkos::TPINode"
    //    NumThreads   = 8
    // 
    // In this scenario, nodes that are equivalent to zero module 2, i.e., even nodes, will be selected to use ThrustGPUNode objects
    // and initialized with the parameter list containing
    //    NodeType   = "Kokkos::ThrustGPUNode"
    //    DeviceNumber = 0
    //    Verbose      = 1
    // Nodes that are equivalent to one modulo 2, i.e., odd nodes, will be selected to use TPINode objects and initialized with the 
    // parameter list containing 
    //    NodeType   = "Kokkos::TPINode"
    //    NumThreads = 8
    // 
    // If multiple node designation sublists match the processor rank, then the first enounctered node designation will be used.
    // I don't know if ParameterList respects any ordering, therefore, multiple matching designations are to be avoided.

    const int myrank = comm_->getRank();
    std::string desigNode("");
    bool matchFound = false;
    for (Teuchos::ParameterList::ConstIterator it = pl.begin(); it != pl.end(); ++it) {
      if (it->second.isList()) {
        int parsedLen, M, N;
        const std::string &name = it->first;
        const Teuchos::ParameterList &sublist = Teuchos::getValue<Teuchos::ParameterList>(it->second);
        // select and assign instList_;
        parsedLen = 0;
        if (std::sscanf(name.c_str(),"%%%d=%d%n",&M,&N,&parsedLen) == 2 && (size_t)parsedLen == name.length()) {
          if ((myrank % M) == N) {
            matchFound = true;
          }
        }
        parsedLen = 0;
        if (std::sscanf(name.c_str(),"=%d%n",&N,&parsedLen) == 1 && (size_t)parsedLen == name.length()) {
          if (myrank == N) {
            matchFound = true;
          }
        }
        parsedLen = 0;
        if (std::sscanf(name.c_str(),"[%d,%d]%n",&M,&N,&parsedLen) == 2 && (size_t)parsedLen == name.length()) {
          if (M <= myrank && myrank <= N) {
            matchFound = true;
          }
        }
        if (name == "default") {
          matchFound = true;
        }
        if (matchFound) {
          try {
            desigNode = sublist.get<std::string>("NodeType");
          }
          catch (Teuchos::Exceptions::InvalidParameterName &e) {
            TEST_FOR_EXCEPTION_PURE_MSG(true, std::runtime_error, 
              std::endl << Teuchos::typeName(*this) << ": Invalid machine file." << std::endl 
              << "Missing parameter \"NodeType\" on Node " << myrank << " for Node designator " << "\"" << name << "\":" << std::endl 
              << sublist << std::endl);
          }
          if (desigNode == "Kokkos::SerialNode") {
            nodeType_ = SERIALNODE;
          }
#ifdef HAVE_KOKKOS_THREADPOOL
          else if (desigNode == "Kokkos::TPINode") {
            nodeType_ = TPINODE;
          }
#endif
#ifdef HAVE_KOKKOS_TBB
          else if (desigNode == "Kokkos::TBBNode") {
            nodeType_ = TBBNODE;
          }
#endif
#ifdef HAVE_KOKKOS_THRUST
          else if (desigNode == "Kokkos::ThrustGPUNode") {
            nodeType_ = THRUSTGPUNODE;
          }
#endif
          else {
            matchFound = false;
          }
          if (matchFound) {
            instList_ = sublist;
            break;
          }
        }
      }
    }
    if (!matchFound) {
      TEST_FOR_EXCEPTION_PURE_MSG(true, std::runtime_error, 
          Teuchos::typeName(*this) << ": No matching node type on rank " << myrank);
    }
  } 

  HybridPlatform::~HybridPlatform() 
  {}

  Teuchos::RCP<const Teuchos::Comm<int> > 
  HybridPlatform::getComm() const {
    return comm_;
  }

  void HybridPlatform::createNode() {
    using Teuchos::rcp;
    if (nodeCreated_) return;
    switch (nodeType_) {
      case SERIALNODE:
        serialNode_ = rcp(new Kokkos::SerialNode(instList_));
        break;
#ifdef HAVE_KOKKOS_TBB
      case TBBNODE:
        tbbNode_ = rcp(new Kokkos::TBBNode(instList_));
        break;
#endif        
#ifdef HAVE_KOKKOS_THREADPOOL
      case TPINODE:
        tpiNode_  = rcp(new Kokkos::TPINode(instList_));
        break;
#endif        
#ifdef HAVE_KOKKOS_THRUST
      case THRUSTGPUNODE:
        thrustNode_ = rcp(new Kokkos::ThrustGPUNode(instList_));
        break;
#endif        
      default:
        TEST_FOR_EXCEPTION(true, std::runtime_error, 
            Teuchos::typeName(*this) << "::runUserCode(): Invalid node type." << std::endl);
    } // end of switch
    nodeCreated_ = true;
  }

  template<template<class Node> class UserCode>
  void HybridPlatform::runUserCode() {
    createNode();
    switch (nodeType_) {
      case SERIALNODE:
        UserCode<Kokkos::SerialNode>::run(instList_,comm_, serialNode_);
        break;
#ifdef HAVE_KOKKOS_TBB
      case TBBNODE:
        UserCode<Kokkos::TBBNode>::run(instList_,comm_, tbbNode_);
        break;
#endif        
#ifdef HAVE_KOKKOS_THREADPOOL
      case TPINODE:
        UserCode<Kokkos::TPINode>::run(instList_,comm_, tpiNode_);
        break;
#endif        
#ifdef HAVE_KOKKOS_THRUST
      case THRUSTGPUNODE:
        UserCode<Kokkos::ThrustGPUNode>::run(instList_,comm_, thrustNode_);
        break;
#endif        
      default:
        TEST_FOR_EXCEPTION(true, std::runtime_error, 
            Teuchos::typeName(*this) << "::runUserCode(): Invalid node type." << std::endl);
    } // end of switch
  }

} // namespace Tpetra

#endif // TPETRA_HYBRIDPLATFORM_HPP
