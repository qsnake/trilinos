// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceDefs.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFSimpleAddedOpDecl.hpp"
#include "TSFSimpleComposedOpDecl.hpp"
#include "TSFSimpleBlockOpDecl.hpp"
#include "TSFSimpleScaledOpDecl.hpp"
#include "TSFSimpleZeroOpDecl.hpp"
#include "TSFSimpleIdentityOpDecl.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "TSFPreconditionerFactory.hpp"
#include "TSFPreconditioner.hpp"
#include "TSFGenericLeftPreconditioner.hpp"
#include "TSFGenericRightPreconditioner.hpp"
#include "TSFProductVectorSpaceDecl.hpp"
#include "TSFGMRESSolver.hpp"
#include "TSFBICGSTABSolverDecl.hpp"
#include "TSFNOXSolver.H"
#include "TSFLinearSolverBuilder.hpp"
#include "TSFLinearCombinationDecl.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PyTeuchos_Utils.hpp"
//#include "PySundanceNOXSolverHandle.hpp"
#include "PySundanceLinearSolver.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorSpaceImpl.hpp"
#include "TSFVectorImpl.hpp"
#include "TSFBICGSTABSolverImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#include "TSFProductVectorSpaceImpl.hpp"
#include "TSFSimpleBlockOpImpl.hpp"
#include "TSFSimpleZeroOpImpl.hpp"
#include "TSFSimpleIdentityOpImpl.hpp"
#include "TSFSimpleAddedOpImpl.hpp"
#include "TSFSimpleComposedOpImpl.hpp"
#include "TSFSimpleScaledOpImpl.hpp"
#endif

#include "TSFHack.hpp"

  %}




// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


%template(doubleVector) std::vector<double>;

%rename(Vector) Vec;
%rename(VectorType) VecType;
%rename(VectorSpace) VecSpace;
%rename(LinearOperator) LinOp;
%rename(LinearSolver) LinSol;
%rename(NonlinearOperator) NonlinOp;
//%rename(NOXSolver) makeNOXSolver;
%rename(Preconditioner) Precond;
%rename(PreconditionerFactory) PrecondFactory;



/* --------- vector space ------------ */
namespace TSFExtended
{
  template <class Scalar> class Vector;
  template <class Scalar>
  class VectorSpace
  {
  public:
    Vector<Scalar> createMember();

    int dim() const ;
  };

  %template(VecSpace) VectorSpace<double>;
}

/* --------- vector ------------ */
namespace TSFExtended
{
  template <class Scalar> class Vector
  {
  public:
    Vector();
    ~Vector();

    void setElement(int globalIndex, const Scalar& x); 

    VectorSpace<Scalar> space() const ;

    Vector<Scalar> copy() const ;

    Vector<Scalar> acceptCopyOf(const Vector<Scalar>& x);

    Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

    Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

    Vector<Scalar> reciprocal() const ;

    Vector<Scalar> abs() const ;

    void setToConstant(const Scalar& alpha) ;

    Scalar norm1() const ;

    Scalar norm2() const ;

    Scalar normInf() const ;

    void zero();

    Scalar max() const;

    Scalar max(int& index)const;

    Scalar max(const Scalar& bound, int& index) const;

    Scalar min()const;

    Scalar min(int& index)const;

    Scalar min(const Scalar& bound, int& index) const;

    Vector<Scalar> getBlock(int i) const  ;

    void setBlock(int i, const Vector<Scalar>& x) ;


    %extend 
    {
      int numBlocks() const
      {
        return self->space().numBlocks();
      }

      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }

      Vector<Scalar> __add__(const Vector<Scalar>& other) 
      {
        return (*self) + other;
      }

      Vector<Scalar> __sub__(const Vector<Scalar>& other) 
      {
        return (*self) - other;
      }

      Vector<Scalar> __mul__(const Scalar& other) 
      {
        return (*self) * other;
      }

      Vector<Scalar> __div__(const Scalar& other) 
      {
        return (*self) *(1.0/ other);
      }

      Vector<Scalar> __rmul__(const Scalar& other) 
      {
        return (*self) * other;
      }

      Scalar __mul__(const Vector<Scalar>& other) 
      {
        return (*self) * other;
      }

      
      Scalar __getitem__(int globalIndex) const 
      {
        return self->getElement(globalIndex);
      }
      
      void __setitem__(int globalIndex, const Scalar& value)
      {
        return self->setElement(globalIndex, value);
      }
    }
  };

  %template(Vec) Vector<double>;

}

/* --------- vector space ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class VectorSpace
  {
  public:
    VectorSpace();
    ~VectorSpace();

    Vector<Scalar> createMember();


    /** return the number of subblocks. */
    int numBlocks() const ;

    /** get the i-th subblock */
    VectorSpace<Scalar> getBlock(int i) const ;


    /** set the i-th subblock */
    void setBlock(int i, const VectorSpace<Scalar>& space);

    %extend 
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(VecSpace) VectorSpace<double>;
}


%rename(BlockVectorSpace) makeBlockVectorSpace;

%inline %{
  /* Create block vector space */
  TSFExtended::VectorSpace<double> 
    makeBlockVectorSpace(TSFExtended::VectorSpace<double> vs)
  {
    return TSFExtended::productSpace(vs);
  }

  /* Create block vector space */
  TSFExtended::VectorSpace<double> 
    makeBlockVectorSpace(TSFExtended::VectorSpace<double> vs1,
                         TSFExtended::VectorSpace<double> vs2)
  {
    return TSFExtended::productSpace(vs1, vs2);
  }

  /* Create block vector space */
  TSFExtended::VectorSpace<double> 
    makeBlockVectorSpace(TSFExtended::VectorSpace<double> vs1,
                         TSFExtended::VectorSpace<double> vs2,
                         TSFExtended::VectorSpace<double> vs3)
  {
    return TSFExtended::productSpace(vs1, vs2, vs3);
  }

  %}


/* --------- vector type ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class VectorType
  {
  public:
    ~VectorType();
    VectorType();

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }

    %extend
    {
      VectorSpace<Scalar> createEvenlyPartitionedSpace(int nLocal) const 
      {
        return self->createEvenlyPartitionedSpace(MPIComm::world(), nLocal);
      }
    }
  };

  %template(VecType) VectorType<double>;

}



/* --------- vector type ------------ */
namespace TSFExtended
{
  enum SolverStatusCode {SolveCrashed, SolveFailedToConverge, SolveConverged};
  
  template <class Scalar>
  class SolverState
  {
  public:
    SolverState();
    SolverState(SolverStatusCode finalState, const std::string& msg, 
                int finalIters, const Scalar& finalResid);
    ~SolverState();
    
    std::string stateDescription() const ;

    /** */
    const Scalar& finalResid() const ;

    /** */
    int finalIters() const ;

    /** */
    const SolverStatusCode& finalState() const ;

    /** */
    const std::string& finalMsg() const ;

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        os << *self;
        rtn = os.str();
        return rtn;
      }
    }
  };

}



%rename(EpetraVectorType) makeEpetraVectorType;

%inline %{
  /* Create an epetra vector type */
  TSFExtended::VectorType<double> makeEpetraVectorType()
  {
    return TSFExtended::VectorType<double>(new TSFExtended::EpetraVectorType());
  }
  %}


/* --------- linear operator ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class LinearOperator
  {
  public:
    LinearOperator();
    ~LinearOperator();

    /** Return the domain */
    const VectorSpace<Scalar> domain() const ;

    /** Return the range */
    const VectorSpace<Scalar> range() const ;


    /** return number of block rows */
    int numBlockRows() const;
      

    /** return number of block cols */
    int numBlockCols() const;
      

    /** get the (i,j)-th block */
    LinearOperator<Scalar> getBlock(const int &i, const int &j) const ;

    /** set the (i,j)-th block 
     *  If the domain and/or the range are not set, then we
     *  are building the operator
     */
    void setBlock(int i, int j, 
                  const LinearOperator<Scalar>& sub);

    /**
     * Return a TransposeOperator.
     */
    LinearOperator<Scalar> transpose() const ; 


    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }

      Vector<Scalar> __mul__(const Vector<Scalar>& in) const 
      {
        Vector<Scalar> out;
        self->apply(in, out);
        return out;
      }

      LinearOperator<Scalar> __mul__(const LinearOperator<Scalar>& other) const 
      {
        return (*self) * other;
      }

      LinearOperator<Scalar> __add__(const LinearOperator<Scalar>& other) const 
      {
        return (*self) + other;
      }

      LinearOperator<Scalar> __sub__(const LinearOperator<Scalar>& other) const 
      {
        return (*self) + -1.0*other;
      }
      
      LinearOperator<Scalar> __mul__(const Scalar& other) const 
      {
        return other * (*self);
      }
      
      LinearOperator<Scalar> __rmul__(const Scalar& other) const 
      {
        return other * (*self);
      }

      
    }
  };

  %template(LinOp) LinearOperator<double>;

}

%inline %{
  /* Create block operator */
  TSFExtended::LinearOperator<double> 
    BlockOperator(const TSFExtended::VectorSpace<double>& domain,
      const TSFExtended::VectorSpace<double>& range)
  {
    return makeBlockOperator(domain, range);
  }

  %}



%rename(IdentityOperator) makeIdentityOperator;

%inline %{
  /* Create block operator */
  TSFExtended::LinearOperator<double> 
    makeIdentityOperator(const TSFExtended::VectorSpace<double>& space)
  {
    return identityOperator(space);
  }

  %}

%rename(ZeroOperator) makeZeroOperator;

%inline %{
  /* Create block operator */
  TSFExtended::LinearOperator<double> 
    makeZeroOperator(const TSFExtended::VectorSpace<double>& domain,
      const TSFExtended::VectorSpace<double>& range)
  {
    return zeroOperator(domain, range);
  }

  %}


/* --------- linear solver ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class LinearSolver
  {
  public:
    LinearSolver();
    ~LinearSolver();

    SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                              const Vector<Scalar>& rhs,
                              Vector<Scalar>& soln) const ;

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(LinSol) LinearSolver<double>;


 //  %extend LinOp 
//   {
//     LinOp inverse(const LinearSolver<Scalar>& solver) const
//     {
//       return LinOp(new InverseOperator<double>(*self, solver));
//     }
//   }

}



%rename(BICGSTABSolver) makeBICGSTABSolver;
%rename(GMRESSolver) makeGMRESSolver;

%inline %{
  TSFExtended::LinearSolver<double> makeGMRESSolver(const Teuchos::ParameterList& params)
  {
    return LinearSolver<double>(new TSFExtended::GMRESSolver<double>(params));
  }
  %}

%inline %{
  TSFExtended::LinearSolver<double> makeBICGSTABSolver(const Teuchos::ParameterList& params)
  {
    return LinearSolver<double>(new TSFExtended::BICGSTABSolver<double>(params));
  }
  %}

%inline %{
  TSFExtended::LinearSolver<double> makeGMRESSolver(const Teuchos::ParameterList& params,
                                                    const TSFExtended::PreconditionerFactory<double>& precond)
  {
    return LinearSolver<double>(new TSFExtended::GMRESSolver<double>(params, precond));
  }
  %}

%inline %{
  TSFExtended::LinearSolver<double> makeBICGSTABSolver(const Teuchos::ParameterList& params,
                                                    const TSFExtended::PreconditionerFactory<double>& precond)
  {
    return LinearSolver<double>(new TSFExtended::BICGSTABSolver<double>(params, precond));
  }
  %}


 

%inline %{
  /* Read a linear solver from an XML file */
  TSFExtended::LinearSolver<double> readSolver(const std::string& filename)
  {
    Teuchos::ParameterXMLFileReader reader(filename);
    Teuchos::ParameterList solverParams = reader.getParameters();
    TSFExtended::LinearSolver<double> solver 
      = TSFExtended::LinearSolverBuilder::createSolver(solverParams);
    return solver;
  }
  %}

%inline %{
  /* Read a linear solver from a parameter list */
  TSFExtended::LinearSolver<double> buildSolver(const Teuchos::ParameterList& params)
  {
    TSFExtended::LinearSolver<double> solver ;
    try
      {
        solver = TSFExtended::LinearSolverBuilder::createSolver(params);
      }
    catch(std::exception& e)
      {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        std::cerr << "detected exception "
                  << e.what() << " in buildSolver()" << std::endl;
      }
    return solver;
  }
  %}


/* --------- preconditioner ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class Preconditioner
  {
  public:
    Preconditioner();
    ~Preconditioner();

    
    
    /** Left preconditioner */
    LinearOperator<Scalar> left() const ;
    
    /** Right preconditioner */
    LinearOperator<Scalar> right() const ;
    
    /** return true if this preconditioner has both left and
     * right components. */
    bool isTwoSided() const ;
    
    /** return true if this preconditioner has a nontrivial left component */
    bool hasLeft() const ;
    
    /** return true if this preconditioner has
     * a nontrivial right component */
    bool hasRight() const ;
    
    /** return true if this preconditioner has neither left nor
     * right operators defined */
    bool isIdentity() const ;

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(Precond) Preconditioner<double>;
}


/* --------- preconditioner factory ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class PreconditionerFactory
  {
  public:
    PreconditionerFactory();
    ~PreconditionerFactory();

    Preconditioner<Scalar> createPreconditioner(const LinearOperator<Scalar>& A) const ;

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(PrecondFactory) PreconditionerFactory<double>;
}


%rename(ILUKPreconditionerFactory) makeILUKPreconditionerFactory;
%inline %{
  /* Create ILUK preconditioner factory  */
  TSFExtended::PreconditionerFactory<double> 
    makeILUKPreconditionerFactory(const Teuchos::ParameterList& params)
  {
    return TSFExtended::PreconditionerFactory<double>(new TSFExtended::ILUKPreconditionerFactory<double>(params));
  }
  %}

%rename(GenericLeftPreconditioner) makeGenericLeftPreconditioner;

%inline %{
  /* Create generic left preconditioner  */
  TSFExtended::Preconditioner<double> 
    makeGenericLeftPreconditioner(const TSFExtended::LinearOperator<double>& left)
  {
    return TSFExtended::Preconditioner<double>(new TSFExtended::GenericLeftPreconditioner<double>(left));
  }
  %}

%rename(GenericRightPreconditioner) makeGenericRightPreconditioner;

%inline %{
  /* Create generic left preconditioner  */
  TSFExtended::Preconditioner<double> 
    makeGenericRightPreconditioner(const TSFExtended::LinearOperator<double>& op)
  {
    return TSFExtended::Preconditioner<double>(new TSFExtended::GenericRightPreconditioner<double>(op));
  }
  %}





/* --------- nonlinear operator ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class NonlinearOperator
  {
  public:
    NonlinearOperator();
    ~NonlinearOperator();

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(NonlinOp) NonlinearOperator<double>;

}

namespace NOX
{
  namespace StatusTest
  {
    enum StatusType {Unevaluated, Unconverged, Converged, Failed};
  }
}


namespace TSFExtended
{
class NOXSolver 
{
public:
  /** */
  NOXSolver();
  /** */
  NOXSolver(const Teuchos::ParameterList& params);

  %extend {
    NOXSolver(PyObject* dict)
    {
      Teuchos::ParameterList params = dict2ParameterList(dict);
      return new NOXSolver(params);
    }
    }
  
  /** */
  NOX::StatusTest::StatusType solve(const NonlinearOperator<double>& F, 
    Vector<double>& soln) const ;

  /** */
  const LinearSolver<double>& linSolver() const ;

};

%extend NOXSolver {
  NOXSolver(PyObject* dict)
  {
    Teuchos::ParameterList params = dict2ParameterList(dict);
    return NOXSolver(params);
  }
}
}



%inline %{
  /* Create a nonlinear solver from a Python dictionary */
  TSFExtended::NOXSolver makeNOXSolver(PyObject* dict)
  {
    Teuchos::ParameterList params = dict2ParameterList(dict);
    return NOXSolver(params);
  }
  /* Create a nonlinear solver from a Python dictionary */
  TSFExtended::NOXSolver makeNOXSolver(const Teuchos::ParameterList& params)
  {
    return NOXSolver(params);
  }
  %}




%inline %{
  namespace TSFExtended
  {
    SolverState<double> 
    PySundanceLinearSolver_solve(const PySundanceLinearSolver* solver,
                                 const LinearOperator<double>& op,
                                 const Vector<double>& rhs,
                                 Vector<double>& soln)
    {
      swig_type_info* opType = SWIG_TypeQuery("TSFExtended::LinearOperator<double>*");
      TEST_FOR_EXCEPTION(opType==0, runtime_error,
                         "swig could not find a match for type name "
                         "[TSFExtended::LinearOperator<double>]");


      swig_type_info* vecType = SWIG_TypeQuery("TSFExtended::Vector<double>*");
      TEST_FOR_EXCEPTION(vecType==0, runtime_error,
                         "swig could not find a match for type name "
                         "[TSFExtended::Vector<double>]");


      swig_type_info* stateType = SWIG_TypeQuery("TSFExtended::SolverState<double>*");
      TEST_FOR_EXCEPTION(stateType==0, runtime_error,
                         "swig could not find a match for type name "
                         "[TSFExtended::SolverState<double>]");


      PyObject* opObj = SWIG_NewPointerObj( (void*) &op, opType, 0);
      PyObject* rhsObj = SWIG_NewPointerObj( (void*) &rhs, vecType, 0);
      PyObject* x0Obj = SWIG_NewPointerObj( (void*) &soln, vecType, 0);

      PyObject* result = solver->pySolve(opObj, rhsObj, x0Obj);

      if (0 == result) {
        PyErr_Print();
        return SolverState<double>(SolveCrashed, "null result from PySundanceLinearSolver",
                                   1, 0.0);
      }

      PyObject* solnObj = 0;
      PyObject* stateObj = 0 ;

      Vector<double>* x = 0 ;
      SolverState<double>* state = 0 ;

      int isTuple = PyTuple_Check(result);

      if (isTuple)
        {
          int size = PyTuple_Size(result);
          switch(size)
            {
            case 2:
              stateObj = PyTuple_GetItem(result, 1);
              TEST_FOR_EXCEPTION(stateObj==0, runtime_error,
                                 "null solver state in PySundanceLinearSolver_solve()");
              SWIG_Python_ConvertPtr(stateObj, (void**) &state, stateType,  
                                     SWIG_POINTER_EXCEPTION | 0);
            case 1:
              solnObj = PyTuple_GetItem(result, 0);
              TEST_FOR_EXCEPTION(solnObj==0, runtime_error,
                                 "null solution object in PySundanceLinearSolver_solve()");
              SWIG_Python_ConvertPtr(solnObj, (void**) &x, vecType,  
                                     SWIG_POINTER_EXCEPTION | 0);
              break;
            default:
              TEST_FOR_EXCEPTION(size < 1 || size > 2, runtime_error,
                                 "invalid return value size " << size 
                                 << " in PySundanceLinearSolver_solve()");
            }
        }
      else
        {
          SWIG_Python_ConvertPtr(result, (void**) &x, vecType,  
                                 SWIG_POINTER_EXCEPTION | 0);
        }

      TEST_FOR_EXCEPTION(x==0, runtime_error, "null return vector in "
                         " PySundanceLinearSolver_solve()");
      soln = *x;

      SolverState<double> rtn(SolveConverged, "unknown solve state", 1, 0);
      if (state!=0) rtn = *state;
      

      Py_DECREF(result); // All done with returned result object

      return rtn;
    }
  }
  %}



