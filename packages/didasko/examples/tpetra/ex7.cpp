// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
//
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_ScalarTraits.hpp"

// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 28-Nov-05

const int MyScalarTypeSize = 2;
const int MyScalarTypePackSize = MyScalarTypeSize * sizeof(double) / sizeof(char);

class MyScalarType
{
public:
  inline MyScalarType(const double value = 0.0)
  {
    for (int i = 0 ; i < MyScalarTypeSize ; ++i)
      data[i] = value;
  }

  inline MyScalarType& operator=(const double rhs)
  {
    for (int i = 0 ; i < MyScalarTypeSize ; ++i)
      data[i] = rhs;
    return *this;
  }

  inline MyScalarType& operator=(const MyScalarType rhs)
  {
    for (int i = 0 ; i < MyScalarTypeSize ; ++i)
      data[i] = rhs.data[i];
    return *this;
  }

  inline MyScalarType& operator+=(const MyScalarType rhs)
  {
    for (int i = 0 ; i < MyScalarTypeSize ; ++i)
      data[i] += rhs.data[i];
    return *this;
  }

  inline MyScalarType operator*(const int i) const
  {
    return(MyScalarType(i));
  }

  inline double Magnitude() const
  {
    double res = 0.0;
    for (int i = 0 ; i < MyScalarTypeSize ; ++i)
      res += data[i] * data[i];
    return(sqrt(res));
  }

  double data[MyScalarTypeSize];
};

inline ostream& operator<<(ostream& os, const MyScalarType& obj)
{
  for (int i = 0 ; i < MyScalarTypeSize ; ++i)
    os << obj.data[i] << " ";
  return(os);
}

namespace Teuchos {
  template<>
  struct ScalarTraits<MyScalarType>
  {
    //! Madatory typedef for result of magnitude
    typedef double magnitudeType;

    //! Determines if scalar type is complex
    static const bool isComplex = false;
    //! Determines if scalar type supports relational operators such as <, >, <=, >=.
    static const bool isComparable = false;
    //! Determines if scalar type have machine-specific parameters (i.e. eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax() are supported)
    static const bool hasMachineParameters = false;
    //! Returns relative machine precision.
    static inline magnitudeType eps()   
    { 
      return(Teuchos::ScalarTraits<double>::eps());
    }

    //! Returns safe minimum (sfmin), such that 1/sfmin does not overflow.
    static inline magnitudeType sfmin() 
    { 
      return(Teuchos::ScalarTraits<double>::sfmin());
    }

    //! Returns the base of the machine.
    static inline magnitudeType base()  
    { 
      return(Teuchos::ScalarTraits<double>::base());
    }

    //! Returns \c eps*base.
    static inline magnitudeType prec()  
    { 
      return(Teuchos::ScalarTraits<double>::prec());
    }

    //! Returns the number of (base) digits in the mantissa.
    static inline magnitudeType t()     
    { 
      return(Teuchos::ScalarTraits<double>::t());
    }

    //! Returns 1.0 when rounding occurs in addition, 0.0 otherwise
    static inline magnitudeType rnd()   
    { 
      return(Teuchos::ScalarTraits<double>::rnd());
    }

    //! Returns the minimum exponent before (gradual) underflow.
    static inline magnitudeType emin()  
    { 
      return(Teuchos::ScalarTraits<double>::emin());
    }

    //! Returns the underflow threshold - \c base^(emin-1)
    static inline magnitudeType rmin()  
    { 
      return(Teuchos::ScalarTraits<double>::rmin());
    }

    //! Returns the largest exponent before overflow.
    static inline magnitudeType emax()  
    { 
      return(Teuchos::ScalarTraits<double>::emin());
    }

    //! Overflow theshold - \c (base^emax)*(1-eps)
    static inline magnitudeType rmax()  
    { 
      return(Teuchos::ScalarTraits<double>::rmax());
    }

    //! Returns the magnitudeType of the scalar type \c a.
    static inline magnitudeType magnitude(MyScalarType a) 
    { 
      return(a.Magnitude());
    }

    //! Returns representation of zero for this scalar type.
    static inline MyScalarType zero()                     
    { 
      return(MyScalarType(0.0));
    }

    //! Returns representation of one for this scalar type.
    static inline MyScalarType one()                
    { 
      return(MyScalarType(1.0));
    }

    //! Returns the conjugate of the scalar type \c a.
    static inline MyScalarType conjugate(MyScalarType a) 
    {
      return(a);
    }

    //! Returns a number that represents NaN.
    static inline MyScalarType nan()                      
    { 
      return(Teuchos::ScalarTraits<double>::nan());
    }

    //! Returns <tt>true</tt> if <tt>x</tt> is NaN or Inf.
    static inline bool isnaninf(const MyScalarType& x)     
    { 
      return(false);
    }

    //! Seed the random number generator returned by <tt>random()</tt>.
    static inline void seedrandom(unsigned int s) 
    {
      Teuchos::ScalarTraits<double>::seedrandom(s);
    }
  
    //! Returns a random number (between -one() and +one()) of this scalar type.
    static inline MyScalarType random()                   
    {
      return(MyScalarType(Teuchos::ScalarTraits<double>::random()));
    }

    //! Returns the name of this scalar type.
    static inline std::string name()           
    {
      return("MyScalarType");
    }

    //! Returns a number of magnitudeType that is the square root of this scalar type \c x. 
    static inline magnitudeType squareroot(MyScalarType x) 
    { 
      throw(1);
    }
  };

} // namespace Teuchos


#ifdef TPETRA_MPI

namespace Tpetra 
{
  template<>
  struct MpiTraits<MyScalarType> 
  {
    static inline MPI_Datatype datatype() {return(MPI_CHAR);};
    static inline int count(int userCount) 
    {
      return(userCount * (MyScalarTypePackSize));
    }

    static MPI_Op sumOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)MyScalarTypeSumOp, true, &myOp);
      return(myOp);
    };
    static MPI_Op maxOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)MyScalarTypeMaxOp, true, &myOp);
      return(myOp);
    };
    static MPI_Op minOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)MyScalarTypeMinOp, true, &myOp);
      return(myOp);
    };

    static void MyScalarTypeSumOp(MyScalarType* in, MyScalarType* inout, 
                                  int* len, MPI_Datatype* dptr) 
    {
      for(int i = 0; i < *len / MyScalarTypePackSize ; i++) 
      {
        inout[i] += in[i];
      }
    };

    static void MyScalarTypeMaxOp(MyScalarType* in, MyScalarType* inout, 
                                  int* len, MPI_Datatype* dptr) 
    {
      for(int i = 0; i < *len / MyScalarTypePackSize ; i++) 
      {
        if (in[i].Magnitude() > inout[i].Magnitude())
          inout[i] = in[i];
      }
    };

    static void MyScalarTypeMinOp(MyScalarType* in, MyScalarType* inout, 
                                  int* len, MPI_Datatype* dptr) 
    {
      for(int i = 0; i < *len / MyScalarTypePackSize ; i++) 
      {
        if (in[i].Magnitude() < inout[i].Magnitude())
          inout[i] = in[i];
      }
    };
  };
} // namespace Tpetra

#endif // TPETRA_MPI

typedef int OrdinalType;
typedef MyScalarType ScalarType;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  if (Comm.getNumImages() != 2)
  {
    cout << "This example can be ran with 2 processors only" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

  // define working arrays

  OrdinalType size = 2;
  vector<ScalarType> V(size), W(size);

  // define all elements on processor 0, then broadcast
  
  if (Comm.getMyImageID() == 0)
  {
    for (OrdinalType i = OrdinalZero ; i < size ; ++i)
      V[i] = ScalarOne * i;
  }

  int RootImage = 0;
  Comm.broadcast(&V[0], size, RootImage);

  for (OrdinalType i = OrdinalZero ; i < size ; ++i)
  {
    cout << "After broadcast(), on image " << Comm.getMyImageID();
    cout << ", V[" << i << "] = " << V[i] << endl;
  }
  
  if (Comm.getMyImageID() == 0)
  {
    cout << "Now setting the elements so that each processor contains its image ID" << endl;
  }

  // sum operations (equivalently, max/min)

  for (OrdinalType i = OrdinalZero ; i < size ; ++i)
    V[i] = Comm.getMyImageID();

  Comm.sumAll(&V[0], &W[0], size);

  for (OrdinalType i = OrdinalZero ; i < size ; ++i)
  {
    cout << "After maxAll(), on image " << Comm.getMyImageID();
    cout << ", W[" << i << "] = " << W[i] << endl;
  }

  // HERE THE FUN COMES... (well, sort of...)
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // now define two maps, defined as follows:
  // - SourceMap: 
  //   *) processor 0 owns elements 0, 1, 2
  //   *) processor 1 owns elements 3, 4, 5
  // - TargetMap:
  //   *) processor 0 owns elements 0, 2, 4
  //   *) processor 1 owns elements 1, 3, 5

  vector<int> MySourceGlobalElements(4);
  if (Comm.getMyImageID() == 0)
  {
    MySourceGlobalElements[0] = 0;
    MySourceGlobalElements[1] = 1;
    MySourceGlobalElements[2] = 2;
    MySourceGlobalElements[2] = 3;
  }
  else
  {
    MySourceGlobalElements[0] = 2;
    MySourceGlobalElements[1] = 3;
    MySourceGlobalElements[2] = 4;
    MySourceGlobalElements[3] = 5;
  }

  Tpetra::ElementSpace<OrdinalType> SourceSpace(-1, 4, MySourceGlobalElements, 0, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorSourceSpace(SourceSpace, platformV);
  Tpetra::Vector<OrdinalType, ScalarType> Source(VectorSourceSpace);

  for (int i = 0 ; i < 4 ; ++i)
    Source[i] = MySourceGlobalElements[i];

  cout << Source;

  vector<int> MyTargetGlobalElements(4);
  if (Comm.getMyImageID() == 0)
  {
    MyTargetGlobalElements[0] = 0;
    MyTargetGlobalElements[1] = 2;
    MyTargetGlobalElements[2] = 4;
    MyTargetGlobalElements[3] = 5;
  }
  else
  {
    MyTargetGlobalElements[0] = 1;
    MyTargetGlobalElements[1] = 3;
    MyTargetGlobalElements[2] = 4;
    MyTargetGlobalElements[3] = 5;
  }

  Tpetra::ElementSpace<OrdinalType> TargetSpace(-1, 4, MyTargetGlobalElements, 0, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorTargetSpace(TargetSpace, platformV);
  Tpetra::Vector<OrdinalType, ScalarType> Target(VectorTargetSpace);

  Tpetra::Import<OrdinalType> Importer(SourceSpace, TargetSpace);

  Target.doImport(Source, Importer, Tpetra::Insert);

  cout << Target;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
