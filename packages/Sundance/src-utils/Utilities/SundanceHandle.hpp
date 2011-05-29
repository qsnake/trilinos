/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// **********************************************************************/
 /* @HEADER@ */

#ifndef SundanceHANDLE_HPP
#define SundanceHANDLE_HPP

#include "SundanceDefs.hpp"
#include "SundanceOut.hpp"
#include "SundancePrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceNamedObject.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TypeNameTraits.hpp"


 /** \def This helper macro defines boilerplate constructors for classes deriving
  * from Handle.
  *
  * If class MyHandle is a handle to a type MyType, simply 
  * put
  * \code
  * HANDLE_CTORS(MyHandle, MyType);
  * \endcode
  * in the class declaration of MyHandle and the macro will create 
  * an empty ctor, a ctor from a smart ptr, and a ctor from a raw pointer. 
  * The macro will also create appropriate doxygen for the handle ctors */
#define HANDLE_CTORS(handle, contents) \
  /** Empty ctor */ \
handle() : Sundance::Handle<contents >() {;} \
  /** Construct a #handle with a raw pointer to a #contents */ \
  handle(Sundance::Handleable<contents >* rawPtr) : Sundance::Handle<contents >(rawPtr) {;} \
  /** Construct a #handle with a smart pointer to a #contents */ \
  handle(const Teuchos::RCP<contents >& smartPtr) : Sundance::Handle<contents >(smartPtr){;}





namespace Sundance
{
using namespace Teuchos;

/** This traits class is used to extract the non-const version of
 * a template argument. The generic case returns the template argument. */
template <class X>
class ConstHandleTraits
{
public:
  typedef X NonconstType;
};


/** Specialization of CHT to types "const X". The nonconst type can
 * be extracted from the template argument. */
template <class X>
class ConstHandleTraits<const X>
{
public:
  typedef X NonconstType;
};


/**
 * Class Sundance::Handle provides a general implementation
 * of the common features of reference-counted handles.
 */
template <class PointerType>
class Handle : public ObjectWithVerbosityBase
{
public:
  /** Empty ctor  */
  Handle() : ptr_() {;}

  /** Construct from a smart pointer */
  Handle(const RCP<PointerType>& _ptr) : ptr_(_ptr) {;}

  /** Construct from a raw pointer to a Handleable.  */
  Handle(Handleable<PointerType>* rawPtr) : ptr_(rawPtr->getRcp()) {;}

  /** Read-only access to the underlying smart pointer. */
  const RCP<PointerType>& ptr() const {return ptr_;}

  /** Read-write access to the underlying smart pointer. */
  RCP<PointerType>& ptr() {return ptr_;}

  /** 
   * Print to a stream using the Printable interface. 
   * If the contents of the handle cannot be 
   * downcasted or crosscasted to a Printable*, an exception
   * will be thrown 
   */
  void print(std::ostream& os) const ;


  /** 
   * Return a short descriptive std::string using the Describable interface.
   * If the contents of the handle cannot be 
   * downcasted or crosscasted to a Describable*, an exception
   * will be thrown. 
   */
  std::string description() const ;

  /** */
  void setName(const std::string& name)
    {
      NamedObject* n = dynamic_cast<NamedObject*>(ptr_.get());
      if (n!=0) n->setName(name);
    }

  /** */
  std::string name() const 
    {
      NamedObject* n = dynamic_cast<NamedObject*>(ptr_.get());
      if (n!=0) return n->name();
      return "AnonymousHandle";
    }

  /** 
   * Return the verbosity setting using the ObjectWithVerbosity
   * interface. If the contents of the handle cannot be downcasted
   * or crosscasted into an ObjectWithVerbosity, a value of
   * zero will be returned.
   */
  int verb() const 
    {
      const ObjectWithVerbosityBase* v 
        = dynamic_cast<const ObjectWithVerbosityBase*>(ptr_.get());
        
      if (v) return v->verb();
      return 0;
    }

  /** 
   * Set the verbosity level of the object using the  ObjectWithVerbosity
   * interface. If the contents of the handle cannot be downcasted
   * or crosscasted into an ObjectWithVerbosity, this call will be
   * ignored and a warning printed.
   */
  void setVerbosity(int x) 
    {
      /* Hack warning: this is a trick to deal with the case where
       * PointerType is const, e.g., someone has written a RCP<const X>. 
       * In such a case the cast to a non-const OWVB would fail. 
       * The ConstHandleTraits business lets us extract the underlying type
       * (e.g., X) with which we can do a const cast. */
      typedef typename ConstHandleTraits<PointerType>::NonconstType NC;
      NC* p = const_cast<NC*>(ptr_.get());
      ObjectWithVerbosityBase* v 
        = dynamic_cast<ObjectWithVerbosityBase*>(p);

      if (v) v->setVerbosity(x);
      else
      {
        Out::os() << "WARNING: cannot set verbosity of object=";
        this->print(Out::os());
        Out::os() << std::endl;
      }
    }

  /** Write a fallback description to be used in objects that are
   * neither named, printable, or describable */
  std::string fallbackDescription() const ;

private:
  RCP<PointerType> ptr_;
};

/* implementation of print() */
template <class PointerType> inline 
void Handle<PointerType>::print(std::ostream& os) const 
{
  const NamedObject* n = dynamic_cast<const NamedObject*>(ptr_.get());
  const Printable* p = dynamic_cast<const Printable*>(ptr_.get());
  const Describable* d = dynamic_cast<const Describable*>(ptr_.get());
  const ObjectWithVerbosityBase* v 
    = dynamic_cast<const ObjectWithVerbosityBase*>(ptr_.get());

  if (v != 0)
  {
    if (v->verb() == 0) 
    {
      if (n) os << n->name();
      else if (d) os << d->description();
      else if (p) p->print(os);
      else os << fallbackDescription();
    }
    else if (v->verb()==1)
    {
      if (d) os << d->description();
      else if (p) p->print(os);
      else os << fallbackDescription();
    }
    else
    {
      if (p) p->print(os);
      else os << fallbackDescription();
    }
  }
  else
  {
    if (p!=0) p->print(os);
    else if (d!=0) os << d->description();
    else if (n!=0) os << n->name();
    else os << fallbackDescription();
  }
}

/* implementation of description() */
template <class PointerType> inline
std::string Handle<PointerType>::description() const 
{
  const Describable* d = dynamic_cast<const Describable*>(ptr_.get());
  const NamedObject* n = dynamic_cast<const NamedObject*>(ptr_.get());
  const ObjectWithVerbosityBase* v 
    = dynamic_cast<const ObjectWithVerbosityBase*>(ptr_.get());
  TeuchosOStringStream oss;

  if (v != 0)
  {
    if (v->verb() == 0) 
    {
      if (n) oss << n->name();
      else if (d) oss << d->description();
      else oss << fallbackDescription();
    }
    else
    {
      if (d) oss << d->description();
      else oss << fallbackDescription();
    }
  }
  else
  {
    if (d!=0) oss << d->description();
    else if (n!=0) oss << n->name();
    else oss << fallbackDescription();
  }
  return oss.str();
}

template <class PointerType> inline
std::string Handle<PointerType>::fallbackDescription() const
{
  typedef typename ConstHandleTraits<PointerType>::NonconstType NC;
  TeuchosOStringStream oss;

  oss << "Handle[" << TypeNameTraits<NC>::name()
      << ", ptr=" << ptr_.get() << "]";
  return oss.str();
}



template <class PointerType> inline
std::ostream& operator<<(std::ostream& os, const Sundance::Handle<PointerType>& h)
{
  h.print(os);
  return os;
}

}

#define STREAM_OUT(handleType) \
                        inline std::ostream& operator<<(std::ostream& os, const handleType& h) \
                        {h.print(os); return os;}



#endif

