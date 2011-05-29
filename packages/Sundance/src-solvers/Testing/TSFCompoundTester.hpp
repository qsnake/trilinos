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


#ifndef TSF_COMPOUNDTESTER_HPP
#define TSF_COMPOUNDTESTER_HPP

#include "TSFLinearOperatorDecl.hpp"
#include "TSFSimpleComposedOpDecl.hpp"
#include "TSFSimpleScaledOpDecl.hpp"
#include "TSFSimpleAddedOpDecl.hpp"
#include "TSFSimpleDiagonalOpDecl.hpp"
#include "TSFTesterBase.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "TSFLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFSimpleComposedOpImpl.hpp"
#include "TSFSimpleScaledOpImpl.hpp"
#include "TSFSimpleAddedOpImpl.hpp"
#include "TSFSimpleDiagonalOpImpl.hpp"
#include "TSFRandomSparseMatrixBuilderImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;

using Thyra::TestSpecifier;

namespace TSFExtended
{

/** */
template <class Scalar>
class CompoundTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  CompoundTester(const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& B,
    const TestSpecifier<Scalar>& sumSpec,
    const TestSpecifier<Scalar>& composedSpec,
    const TestSpecifier<Scalar>& scaledSpec,
    const TestSpecifier<Scalar>& diagSpec);

  /** */
  bool runAllTests() const ;

  /** */
  bool sumTest() const ;

  /** */
  bool composedTest() const ;

  /** */
  bool scaledTest() const ;

  /** */
  bool diagTest() const ;


private:

  LinearOperator<Scalar> A_;

  LinearOperator<Scalar> B_;

  TestSpecifier<Scalar> sumSpec_;

  TestSpecifier<Scalar> composedSpec_;

  TestSpecifier<Scalar> scaledSpec_;

  TestSpecifier<Scalar> diagSpec_;

};

template <class Scalar> 
inline CompoundTester<Scalar>
::CompoundTester(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B,
  const TestSpecifier<Scalar>& sumSpec,
  const TestSpecifier<Scalar>& composedSpec,
  const TestSpecifier<Scalar>& scaledSpec,
  const TestSpecifier<Scalar>& diagSpec)
  : TesterBase<Scalar>(), 
    A_(A),
    B_(B),
    sumSpec_(sumSpec),
    composedSpec_(composedSpec),
    scaledSpec_(scaledSpec),
    diagSpec_(diagSpec)
{;}

template <class Scalar> 
inline bool CompoundTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = sumTest() && pass;
  pass = composedTest() && pass;
  pass = scaledTest() && pass;
  pass = diagTest() && pass;

  return pass;
}

template <class Scalar> 
inline bool CompoundTester<Scalar>
::sumTest() const 
{
  if (sumSpec_.doTest())
  {
    Out::root() << "running operator addition test..." << std::endl;
    LinearOperator<Scalar> sum = A_ + B_;

    Vector<Scalar> x = A_.domain().createMember();
    randomizeVec(x);
    Out::root() << "computing y1 = (A+B)*x..." << std::endl;
    Vector<Scalar> y1 = sum*x;
    Out::root() << "computing y2 = A*x + B*x..." << std::endl;
    Vector<Scalar> y2 = A_*x + B_*x;
    
    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
        
    return checkTest(sumSpec_, err, "operator addition");
  }
  Out::root() << "skipping operator addition test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool CompoundTester<Scalar>
::composedTest() const 
{
  if (composedSpec_.doTest())
  {
    Out::root() << "running operator composition test..." << std::endl;
    LinearOperator<Scalar> composed = A_ * B_;

    Vector<Scalar> x = B_.domain().createMember();
    randomizeVec(x);
    Out::root() << "computing y1 = (A*B)*x..." << std::endl;
    Vector<Scalar> y1 = composed*x;
    Out::root() << "computing y2 = B*x..." << std::endl;
    Vector<Scalar> y2 = B_*x;
    Out::root() << "computing y3 = A*y2..." << std::endl;
    Vector<Scalar> y3 = A_*y2;

    ScalarMag err = (y1 - y3).norm2();

    Out::root() << "|y1-y3| = " << err << std::endl;
    return checkTest(composedSpec_, err, "operator composition");
  }
  Out::root() << "skipping operator composition test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool CompoundTester<Scalar>
::scaledTest() const 
{
  if (scaledSpec_.doTest())
  {
    Out::root() << "running operator scaling test..." << std::endl;
    Scalar alpha = sqrt(2.0);
    LinearOperator<Scalar> scaled = alpha*A_;

    Vector<Scalar> x = A_.domain().createMember();
    randomizeVec(x);
    Out::root() << "computing y1 = (alpha*A)*x..." << std::endl;
    Vector<Scalar> y1 = scaled*x;
    Out::root() << "computing y2 = A*x..." << std::endl;
    Vector<Scalar> y2 = A_*x;
    Out::root() << "computing y3 = alpha*y2..." << std::endl;
    Vector<Scalar> y3 = alpha*y2;

    ScalarMag err = (y1 - y3).norm2();

    Out::root() << "|y1-y3| = " << err << std::endl;
    return checkTest(composedSpec_, err, "operator scaling");
  }
  Out::root() << "skipping operator scaling test..." << std::endl;
  return true;
}

  

template <class Scalar> 
inline bool CompoundTester<Scalar>
::diagTest() const 
{
  if (diagSpec_.doTest())
  {
    Out::root() << "running diagonal operator test..." << std::endl;

    Vector<Scalar> x = A_.domain().createMember();
    randomizeVec(x);

    Vector<Scalar> d = A_.domain().createMember();
    randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);

    Out::root() << "computing y1 = D*x..." << std::endl;
    Vector<Scalar> y1 = D*x;
    Out::root() << "computing y2 = d .* x..." << std::endl;
    Vector<Scalar> y2 = x.dotStar(d);

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    return checkTest(diagSpec_, err, "diagonal operator");
  }
  Out::root() << "skipping diagonal operator test..." << std::endl;
  return true;
}

  
  
}
#endif
