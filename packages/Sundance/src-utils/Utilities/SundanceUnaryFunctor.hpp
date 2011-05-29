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

#ifndef SUNDANCE_UNARYFUNCTOR_H
#define SUNDANCE_UNARYFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceFunctorDomain.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  
  

  /**
   * 
   */
  class UnaryFunctor
  {
  public:
    /** ctor */
    UnaryFunctor(const std::string& name, 
                 const RCP<FunctorDomain>& domain 
                 = rcp(new UnboundedDomain())) 
      : name_(name), h_(fdStep()), domain_(domain) {;}

    /** */
    virtual ~UnaryFunctor(){;}

    /** */
    const std::string& name() const {return name_;}

    /** */
    virtual void eval0(const double* const x, 
                       int nx, 
                       double* f) const = 0 ;
    
    /** */
    virtual void eval1(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx) const ;
    
    /** */
    virtual void eval2(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx) const ;
    
    /** */
    virtual void eval3(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx,
                       double* d3f_dxxx) const ;

    

    /** */
    void evalFDDerivs1(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx) const ;
    /** */
    void evalFDDerivs2(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx) const ;
    /** */
    void evalFDDerivs3(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx,
                       double* d3f_dxxx) const ;

    /** */
    bool testDerivs(const double& x, const double& tol) const ;

    /** */
    bool testInvalidValue(const double& xBad) const ;

    /** */
    bool test(int nx, const double& tol) const ;
    
    /** Specify whether we should test for NAN or INFINITE results. */
    static bool& checkResults() {static bool rtn = false; return rtn;}

    static double& fdStep() {static double rtn = 1.0e-3; return rtn;}

    const RCP<FunctorDomain>& domain() const 
    {return domain_;}
  private:
    std::string name_;

    double h_;

    RCP<FunctorDomain> domain_;
  };
}

/** */
#define SUNDANCE_UNARY_FUNCTOR(opName, functorName, description, domain, \
                               funcDefinition, firstDerivDefinition,    \
                               secondDerivDefinition)                   \
  class functorName : public Sundance::UnaryFunctor                \
  {                                                                     \
  public:                                                               \
    /** ctor for description functor */                                 \
    functorName() : Sundance::UnaryFunctor(#opName, rcp(new domain)) {;} \
      /** virtual dtor */                                               \
      virtual ~functorName(){;}                                         \
      /** Evaluate function at an array of values */                    \
      void eval0(const double* const x, int nx, double* f) const ;      \
      /** Evaluate function and first derivative at an array of values */ \
      void eval1(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df) const ;                                    \
      /** Evaluate function and first two derivatives at an array of values */ \
      void eval2(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df,                                            \
                 double* d2f_dxx) const ;                               \
  };                                                                    \
  inline void functorName::eval0(const double* const x, int nx, double* f) const \
  {                                                                     \
    if (checkResults())                                                 \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
          {                                                             \
            f[i] = funcDefinition;                                      \
            TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i]), RuntimeError, "Non-normal floating point result detected in evaluation of unary functor " << name() << " at argument " << x[i]); \
          }                                                             \
     }                                                                  \
   else                                                                 \
     {                                                                  \
       for (int i=0; i<nx; i++) f[i] = funcDefinition;                  \
     }                                                                  \
}                                                                       \
  inline void functorName::eval1(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df) const                        \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
             TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i]), \
                                RuntimeError,                           \
                                "Non-normal floating point result detected in " \
                                "evaluation of unary functor "          \
                                << name() << " at argument " << x[i]);  \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
           }                                                            \
      }                                                                 \
}                                                                       \
  inline void functorName::eval2(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df,                              \
                               double* d2f) const                       \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
              df[i] = firstDerivDefinition;                             \
              d2f[i] = secondDerivDefinition;                           \
              TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i])||Teuchos::ScalarTraits<double>::isnaninf(d2f[i]), \
                                 RuntimeError,                          \
                                 "Non-normal floating point result detected in " \
                                 "evaluation of unary functor "         \
                                 << name() << " at argument " << x[i] ); \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
               df[i] = firstDerivDefinition;                            \
               d2f[i] = secondDerivDefinition;                          \
           }                                                            \
      }                                                                 \
}



/** */
#define SUNDANCE_UNARY_FUNCTOR3(opName, functorName, description, domain, \
                               funcDefinition, firstDerivDefinition,    \
                                secondDerivDefinition, thirdDerivDefinition) \
  class functorName : public Sundance::UnaryFunctor                \
  {                                                                     \
  public:                                                               \
    /** ctor for description functor */                                 \
    functorName() : Sundance::UnaryFunctor(#opName, rcp(new domain)) {;} \
      /** virtual dtor */                                               \
      virtual ~functorName(){;}                                         \
      /** Evaluate function at an array of values */                    \
      void eval0(const double* const x, int nx, double* f) const ;      \
      /** Evaluate function and first derivative at an array of values */ \
      void eval1(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df) const ;                                    \
      /** Evaluate function and first two derivatives at an array of values */ \
      void eval2(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df,                                            \
                 double* d2f_dxx) const ;                               \
      /** Evaluate function and first thress derivatives at an array of values */ \
      void eval3(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df,                                            \
                 double* d2f_dxx,                                       \
                 double* d3f_dxxx) const ;                                        \
  };                                                                    \
  inline void functorName::eval0(const double* const x, int nx, double* f) const \
  {                                                                     \
    if (checkResults())                                                 \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
          {                                                             \
            f[i] = funcDefinition;                                      \
            TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i]), RuntimeError, "Non-normal floating point result detected in evaluation of unary functor " << name() << " at argument " << x[i]); \
          }                                                             \
     }                                                                  \
   else                                                                 \
     {                                                                  \
       for (int i=0; i<nx; i++) f[i] = funcDefinition;                  \
     }                                                                  \
}                                                                       \
  inline void functorName::eval1(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df) const                        \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
             TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i]), \
                                RuntimeError,                           \
                                "Non-normal floating point result detected in " \
                                "evaluation of unary functor "          \
                                << name() << " at argument " << x[i]);  \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
           }                                                            \
      }                                                                 \
}                                                                       \
  inline void functorName::eval2(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df,                              \
                               double* d2f) const                       \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
              df[i] = firstDerivDefinition;                             \
              d2f[i] = secondDerivDefinition;                           \
              TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i])||Teuchos::ScalarTraits<double>::isnaninf(d2f[i]), \
                                 RuntimeError,                          \
                                 "Non-normal floating point result detected in " \
                                 "evaluation of unary functor "         \
                                 << name() << " at argument " << x[i] ); \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
               df[i] = firstDerivDefinition;                            \
               d2f[i] = secondDerivDefinition;                          \
           }                                                            \
      }                                                                 \
}                                                                     \
  inline void functorName::eval3(const double* const x,                 \
                                 int nx,                                \
                                 double* f,                             \
                                 double* df,                            \
                                 double* d2f,                           \
                                 double* d3f) const                     \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
              df[i] = firstDerivDefinition;                             \
              d2f[i] = secondDerivDefinition;                           \
              d3f[i] = thirdDerivDefinition;                           \
              TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i])||Teuchos::ScalarTraits<double>::isnaninf(d2f[i])||Teuchos::ScalarTraits<double>::isnaninf(d3f[i]), \
                                 RuntimeError,                          \
                                 "Non-normal floating point result detected in " \
                                 "evaluation of unary functor "         \
                                 << name() << " at argument " << x[i] ); \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
               df[i] = firstDerivDefinition;                            \
               d2f[i] = secondDerivDefinition;                          \
               d3f[i] = thirdDerivDefinition;                          \
           }                                                            \
      }                                                                 \
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
