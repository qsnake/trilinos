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


#ifndef TSF_TESTERBASE_HPP
#define TSF_TESTERBASE_HPP

#include "TSFLinearOperatorDecl.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "SundanceOut.hpp"

using namespace TSFExtended;
using namespace Teuchos;

using Thyra::TestSpecifier;

namespace TSFExtended
{

  /** */
  template <class Scalar>
  class TesterBase : public ObjectWithClassVerbosity<TesterBase<Scalar> >
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** */
    TesterBase(){;}

    /** */
    virtual ~TesterBase(){;}

    /** */
    virtual bool runAllTests() const = 0 ;


    /** */
    bool checkTest(const TestSpecifier<Scalar>& spec,
                    const ScalarMag& err, 
                    const std::string& testName) const ;

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;

  };

  template <class Scalar> 
  inline void TesterBase<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Thyra::randomize(Scalar(-ST::one()),Scalar(+ST::one()),x.ptr().ptr());
    
  }

  template <class Scalar> 
  inline bool TesterBase<Scalar>
  ::checkTest(const TestSpecifier<Scalar>& spec,
              const ScalarMag& err, 
              const std::string& testName) const 
  {
    bool rtn = true;
    if (err > spec.errorTol())
      {
        Out::root() << testName << " test FAILED: err=" << err << ", tol = " 
             << spec.errorTol() << std::endl;
        rtn = false;
      }
    else if (err > spec.warningTol())
      {
        Out::root() << "WARNING: " << testName << " test err="
             << err << " could not beat tol = " 
             << spec.warningTol() << std::endl;
      }
    else
      {
        Out::root() << "test " << testName << " PASSED with tol=" << spec.errorTol() << std::endl;
      }
    return rtn;
  }
  
}
#endif
