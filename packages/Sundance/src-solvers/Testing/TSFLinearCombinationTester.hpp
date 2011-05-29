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


#ifndef TSF_LINEARCOMBINATIONTESTER_HPP
#define TSF_LINEARCOMBINATIONTESTER_HPP

#include "TSFLinearOperatorDecl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#include "TSFSimpleComposedOpDecl.hpp"
#include "TSFSimpleScaledOpDecl.hpp"
#include "TSFSimpleAddedOpDecl.hpp"
#include "SundanceOut.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFTesterBase.hpp"
#include "TSFRandomSparseMatrixBuilderDecl.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFRandomSparseMatrixBuilderImpl.hpp"
#include "TSFSimpleComposedOpImpl.hpp"
#include "TSFSimpleScaledOpImpl.hpp"
#include "TSFSimpleAddedOpImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;
using namespace Sundance;

using Thyra::TestSpecifier;


#define TESTER(form1, form2)\
  {\
    Out::os() << "testing " #form1 << std::endl;\
    Vector<Scalar> _val1 = form1;\
    Out::os() << "testing " #form2 << std::endl;\
    Vector<Scalar> _val2 = form2;\
    Out::os() << "done testing... checking error" << std::endl;\
    ScalarMag err = (_val1-_val2).norm2();\
    if (!checkTest(spec_, err, "[" #form1 "] == [" #form2 "]")) pass = false;\
  }
    






namespace TSFExtended
{

/** */
template <class Scalar>
class LinearCombinationTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  LinearCombinationTester(int nLocalRows, double onProcDensity, 
    double offProcDensity,
    const VectorType<Scalar>& vecType,
    const TestSpecifier<Scalar>& spec);

  /** */
  bool runAllTests() const ;


private:

  /** */
  bool nonModifyingOpTests() const ;

  /** */
  bool selfModifyingOpTests() const ;

  TestSpecifier<Scalar> spec_;

  int nLocalRows_;
  double onProcDensity_;
  double offProcDensity_;
  VectorType<Scalar> vecType_;

};

template <class Scalar> 
inline LinearCombinationTester<Scalar>
::LinearCombinationTester(int nLocalRows, double onProcDensity, 
  double offProcDensity,
  const VectorType<Scalar>& vecType,
  const TestSpecifier<Scalar>& spec)
  : TesterBase<Scalar>(),
    spec_(spec),
    nLocalRows_(nLocalRows),
    onProcDensity_(onProcDensity),
    offProcDensity_(offProcDensity),
    vecType_(vecType)
{;}

template <class Scalar> 
inline bool LinearCombinationTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = nonModifyingOpTests() && pass;
  pass = selfModifyingOpTests() && pass;

  return pass;
}

template <class Scalar> 
inline bool LinearCombinationTester<Scalar>
::nonModifyingOpTests() const
{
  bool pass = true;

  RandomSparseMatrixBuilder<double> ABuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> BBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> CBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);


  LinearOperator<double> A = ABuilder.getOp();

    
  LinearOperator<Scalar> B = BBuilder.getOp();
  LinearOperator<Scalar> C = CBuilder.getOp();

  Vector<Scalar> x = A.domain().createMember();
  Vector<Scalar> y = A.domain().createMember();
  Vector<Scalar> z = A.domain().createMember();

  randomizeVec(x);
  randomizeVec(y);
  randomizeVec(z);


  TESTER(x*2.0, 2.0*x);

  TESTER(2.0*(x + y), 2.0*x + 2.0*y);

  TESTER(2.0*(x - y), 2.0*x - 2.0*y);

  TESTER((x + y)*2.0, 2.0*x + 2.0*y);

  TESTER((x - y)*2.0, 2.0*x - 2.0*y);

  TESTER(2.0*(x - y), -2.0*(y - x));

  TESTER(0.25*(2.0*(x + y) - 2.0*(x - y)), y);

  TESTER((2.0*A)*x, 2.0*(A*x));

  TESTER(2.0*(A*x), (A*x)*2.0);

  TESTER(A*(B*x), (A*B)*x);

  TESTER(2.0*A*(B*x), A*(B*(2.0*x)));

  TESTER(3.0*(2.0*A)*x, 6.0*(A*x));

  TESTER(y + A*x, A*x + y);


  TESTER(z + (A*x + B*y), (B*y + A*x) + z);


  TESTER(z - (A*x + B*y), -1.0*((B*y + A*x) - z));

  TESTER(C*z + (A*x + B*y), (B*y + A*x) + C*z);

  TESTER(C*z - (A*x + B*y), -1.0*((B*y + A*x) - C*z));

  TESTER(2.0*z + (A*x + B*y), (B*y + A*x) + 2.0*z);

  TESTER(2.0*z - (A*x + B*y), -1.0*((B*y + A*x) - 2.0*z));

  TESTER(A*x - y, -1.0*(y - A*x));

  TESTER(A*x + B*y, B*y + A*x);

  TESTER(A*x - B*y - A*x + B*y +z, z);

  TESTER(2.0*(A*x + y), 2.0*A*x + 2.0*y);

  TESTER(2.0*(A*x + B*y), A*x + B*y + A*x + B*y);

  TESTER(2.0*(y + A*x), 2.0*y + 2.0*(A*x));

  TESTER(x + 2.0*A*y, x + 2.0*(A*y));

  TESTER(2.0*A*y + x, 2.0*(A*y) + x);

  TESTER(2.0*(A*x + B*y), 2.0*A*x + 2.0*B*y);

  TESTER(2.0*(A*x - B*y), 2.0*A*x - 2.0*B*y);

  TESTER(2.0*(A*x + B*y + z), 2.0*A*x + 2.0*B*y + 2.0*z);

  TESTER(2.0*(A*x + 3.0*B*y), 2.0*A*x + 6.0*B*y);

  TESTER(2.0*(A*x + 3.0*(z + B*y)), 2.0*A*x + 6.0*B*y + 6.0*z);

  TESTER(2.0*(z + A*x + B*y + z), 2.0*A*x + 2.0*B*y + 4.0*z);

  TESTER(2.0*(3.0*(z + A*x) + B*y), 6.0*z + 6.0*A*x + 2.0*B*y);

  TESTER(2.0*(3.0*(z + A*x) + 4.0*(B*y + z)), 6.0*z + 6.0*A*x + 8.0*B*y + 8.0*z);
    
  TESTER((A*x + B*y) + (A*y + B*x), (A + B)*x + (A+B)*y);
  TESTER((A*x + B*y) - (A*y + B*x), A*x - A*y + B*y - B*x);

  TESTER((A*x + B*y) + 2.0*(A*y + B*x), A*(x + 2.0*y) + B*(2.0*x + y));
  TESTER((A*x + B*y) - 2.0*(A*y + B*x), A*(x - 2.0*y) + B*(y - 2.0*x));


  return pass;
}


template <class Scalar> 
inline bool LinearCombinationTester<Scalar>
::selfModifyingOpTests() const
{
  bool pass = true;

  RandomSparseMatrixBuilder<double> ABuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> BBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> CBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);

  LinearOperator<double> A = ABuilder.getOp();
  LinearOperator<double> B = BBuilder.getOp();
  LinearOperator<double> C = CBuilder.getOp();

  Vector<Scalar> x = A.domain().createMember();
  Vector<Scalar> y = A.domain().createMember();
  Vector<Scalar> z = A.domain().createMember();

  randomizeVec(x);
  randomizeVec(y);
  randomizeVec(z);

  Vector<Scalar> a = x.copy();
  Vector<Scalar> b = y.copy();
  Vector<Scalar> c = z.copy();
    

  Out::os() << "starting linear combination tests" << std::endl;

  x = 2.0*A*x;
  Scalar err = (x - 2.0*A*a).norm2();
  if (!checkTest(spec_, err, "x=2.0*A*x")) pass = false;

  a = x.copy();
  x = x + y;
  err = (x - (a + y)).norm2();
  if (!checkTest(spec_, err, "x=x+y")) pass = false;

  a = x.copy();
  x = y + x;
  err = (x - (y + a)).norm2();
  if (!checkTest(spec_, err, "x=y+x")) pass = false;

  a = x.copy();
  x = A*x + B*y;
  err = (x - (A*a + B*y)).norm2();
  if (!checkTest(spec_, err, "x=A*x+B*y")) pass = false;

  a = x.copy();
  x = B*y + A*x;
  err = (x - (B*y + A*a)).norm2();
  if (!checkTest(spec_, err, "x=B*y+A*x")) pass = false;

  a = x.copy();
  x = A*x + (B*y + C*x);
  err = (x - (A*a + (B*y + C*a))).norm2();
  if (!checkTest(spec_, err, "x=A*x + (B*y + C*x)")) pass = false;

  a = x.copy();
  x = (A*x + B*y) + C*x;
  err = (x - ((A*a + B*y) + C*a)).norm2();
  if (!checkTest(spec_, err, "x=(A*x + B*y) + C*x")) pass = false;

  /* test assignment of OpTimesLC into empty and non-empty vectors */
  Vector<Scalar> u;
  u = 2.0*A*B*x;
  err = (u - 2.0*A*B*x).norm2();
  if (!checkTest(spec_, err, "(empty) u=2*A*B*x")) pass = false;

  u = 2.0*A*B*x;
  err = (u - 2.0*A*B*x).norm2();
  if (!checkTest(spec_, err, "(non-empty) u=2*A*B*x")) pass = false;

  /* test assignment of LC2 into empty and non-empty vectors */
  Vector<Scalar> v;
  v = 2.0*x + 3.0*y;
  err = (v - (2.0*x + 3.0*y)).norm2();
  if (!checkTest(spec_, err, "(empty) v=2*x + 3*y")) pass = false;

  v = 2.0*x + 3.0*y;
  err = (v - (2.0*x + 3.0*y)).norm2();
  if (!checkTest(spec_, err, "(non-empty) v=2*x + 3*y")) pass = false;


    


  return pass;
}

  
}
#endif
