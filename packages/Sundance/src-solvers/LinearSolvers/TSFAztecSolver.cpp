#include "TSFAztecSolver.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFEpetraMatrix.hpp"
#include "Ifpack_Preconditioner.h"
#include "Ifpack.h"
#include "EpetraTSFOperator.hpp"
#include "Teuchos_basic_oblackholestream.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#endif

#ifdef HAVE_ML
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_operator.h"
#include "ml_aztec_utils.h"
#include "ml_MultiLevelPreconditioner.h"
using namespace ML_Epetra;
#else
#error blarf
#endif

using namespace TSFExtended;
using namespace Teuchos;


AztecSolver::AztecSolver(const ParameterList& params)
  : LinearSolverBase<double>(params),
		options_(AZ_OPTIONS_SIZE),
		parameters_(AZ_PARAMS_SIZE),
    useML_(false),
    useIfpack_(false),
    useUserPrec_(false),
    aztec_recursive_iterate_(false),
    precParams_(),
    userPrec_(),
    aztec_status(AZ_STATUS_SIZE),
    aztec_proc_config(AZ_PROC_SIZE)
{
  setName("AztecSolver");
  initParamMap();

  /* initialize the options and parameters with Aztec's defaults */
	AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));


  /* Set options according to the parameter list */
  ParameterList::ConstIterator iter;
  for (iter=params.begin(); iter != params.end(); ++iter)
  {
    const std::string& name = params.name(iter);
    const ParameterEntry& entry = params.entry(iter);

    if (entry.isList())
    {
      if (name=="Preconditioner")
      {
        precParams_ = params.sublist("Preconditioner");
        TEST_FOR_EXCEPTION(!precParams_.isParameter("Type"), std::runtime_error,
          "preconditioner type not specified in parameter list "
          << precParams_);
        if (precParams_.get<string>("Type") == "ML")
        {
          useML_ = true;
        }
        else if (precParams_.get<string>("Type") == "Ifpack")
        {
          useIfpack_ = true;
        }
	      else if (precParams_.get<string>("Type") == "User")
        {
          useUserPrec_ = true;
        }
        continue;
      }
    }

    /* Check that the param name appears in the table of Aztec params */
    if (paramMap().find(name) == paramMap().end()) continue;

    /* find the integer ID used by Aztec to identify this parameter */
    int aztecCode = paramMap()[name];

    if (name=="Recursive Iterate")
    {
      if (getValue<int>(entry))
        aztec_recursive_iterate_ = true;
      else
        aztec_recursive_iterate_ = false;
      continue;
    }


    /* We now need to figure out what to do with the value of the
     * parameter. If it is a std::string, then it corresponds to a
     * predefined Aztec option value. If it is an integer, then
     * it is the numerical setting for an Aztec option. If it is
     * a double, then it is the numerical setting for an Aztec
     * parameter. */
    if (entry.isType<string>())
    {
      std::string val = getValue<string>(entry);
      TEST_FOR_EXCEPTION(paramMap().find(val) == paramMap().end(),
        std::runtime_error,
        "Aztec solver ctor: [" << val << "] is not a "
        "valid Aztec option value");
      int optionVal = paramMap()[val];
      options_[aztecCode] = optionVal;
    }
    else if (entry.isType<int>())
    {
      int val = getValue<int>(entry);
      options_[aztecCode] = val;
      if (name=="Verbosity") setVerbosity(val);
    }
    else if (entry.isType<double>())
    {
      double val = getValue<double>(entry);
      parameters_[aztecCode] = val;
    }
  }
}


AztecSolver::AztecSolver(const Teuchos::map<int, int>& aztecOptions,
  const Teuchos::map<int, double>& aztecParameters)
  : LinearSolverBase<double>(ParameterList()),
    options_(AZ_OPTIONS_SIZE),
    parameters_(AZ_PARAMS_SIZE),
    useML_(false),
    useIfpack_(false),
    useUserPrec_(false),
    aztec_recursive_iterate_(false),
    precParams_(),
    userPrec_(),
    aztec_status(AZ_STATUS_SIZE),
    aztec_proc_config(AZ_PROC_SIZE)
{
  setName("AztecSolver");
  if (aztecOptions.find(AZ_recursive_iterate) != aztecOptions.end())
  {
    aztec_recursive_iterate_ = true;
  }

  /* initialize the options and parameters with Aztec's defaults */
  AZ_defaults((int*) &(options_[0]), (double*) &(parameters_[0]));

  /* set user-specified options  */
  map<int, int>::const_iterator opIter;
  for (opIter=aztecOptions.begin(); opIter!=aztecOptions.end(); opIter++)
  {
    int opKey = opIter->first;
    if (opKey==AZ_recursive_iterate) continue;
    int opValue = opIter->second;
    options_[opKey] = opValue;
  }

  /* set user-specified params  */
  map<int, double>::const_iterator parIter;
  for (parIter=aztecParameters.begin(); parIter!=aztecParameters.end();
       parIter++)
  {
    int parKey = parIter->first;
    double parValue = parIter->second;
    parameters_[parKey] = parValue;
  }
}


void AztecSolver::updateTolerance(const double& tol)
{
  parameters_[AZ_tol] = tol;
}

SolverState<double> AztecSolver::solve(const LinearOperator<double>& op, 
  const Vector<double>& rhs, 
  Vector<double>& soln) const
{
  RCP<ostream> out;
  if (verb()==0)
  {
    out = rcp(new oblackholestream());
  }
  else
  {
    out = rcp(&Out::os(), false);
  }
  RCP<MultiLevelPreconditioner> mlPrec;
  RCP<Ifpack_Preconditioner> ifpackPrec;
  RCP<Epetra_Operator> prec;

	TSFExtended::Vector<double> bCopy = rhs.copy();
	TSFExtended::Vector<double> xCopy = rhs.copy();

  if (verb() > 4) 
  {
    *out << "rhs=" << bCopy << std::endl;
  }

  Epetra_Vector* b = EpetraVector::getConcretePtr(bCopy);
  Epetra_Vector* x = EpetraVector::getConcretePtr(xCopy);

	Epetra_CrsMatrix& A = EpetraMatrix::getConcrete(op);

  AztecOO aztec(&A, x, b);

  aztec.SetAllAztecOptions((int*) &(options_[0]));
  aztec.SetAllAztecParams((double*) &(parameters_[0]));
  aztec.SetOutputStream(*out);
  aztec.SetErrorStream(*out);
  
  int maxIters = options_[AZ_max_iter];
  double tol = parameters_[AZ_tol];

  if (useML_)
  {
    std::string precType = precParams_.get<string>("Problem Type");
    ParameterList mlParams;
    ML_Epetra::SetDefaults(precType, mlParams);
    //#ifndef TRILINOS_6
    //      mlParams.setParameters(precParams_.sublist("ML Settings"));
    //#else
    ParameterList::ConstIterator iter;
    ParameterList mlSettings = precParams_.sublist("ML Settings");
    for (iter=mlSettings.begin(); iter!=mlSettings.end(); ++iter)
    {
      const std::string& name = mlSettings.name(iter);
      const ParameterEntry& entry = mlSettings.entry(iter);
      mlParams.setEntry(name, entry);
    }
    //#endif
    mlPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(A, mlParams));
    prec = rcp_dynamic_cast<Epetra_Operator>(mlPrec);
  }
  else if (useIfpack_)
  {
    Ifpack precFactory;
    int overlap = precParams_.get<int>("Overlap");
    std::string precType = precParams_.get<string>("Prec Type");

    ParameterList ifpackParams = precParams_.sublist("Ifpack Settings");

    ifpackPrec = rcp(precFactory.Create(precType, &A, overlap));
    prec = rcp_dynamic_cast<Epetra_Operator>(ifpackPrec);
    ifpackPrec->SetParameters(ifpackParams);
    ifpackPrec->Initialize();
    ifpackPrec->Compute();
  }
  else if (useUserPrec_)
  {
    TEST_FOR_EXCEPT(userPrec_.get() == 0);
    prec = userPrec_;
  }
  
  
  if (prec.get() != 0) aztec.SetPrecOperator(prec.get());  
  
  aztec.CheckInput();
  
  /* VEH/RST Parameter to check if we are calling aztec recursively.
   * If so, need to set parameter aztec_recursive_iterate to true. */
  if (aztec_recursive_iterate_)
  {
    aztec.recursiveIterate(maxIters, tol);
  }
  else
  {
    aztec.Iterate(maxIters, tol);
  }
  
  soln = xCopy;

  const double* status = aztec.GetAztecStatus();
  SolverStatusCode state = SolveCrashed;

  std::string msg;
  switch((int) status[AZ_why])
  {
    case AZ_normal:
      state = SolveConverged;
      msg = "converged";
      break;
    case AZ_param:
      state = SolveCrashed;
      msg = "failed: parameter not available";
      break;
    case AZ_breakdown:
      state = SolveCrashed;
      msg = "failed: numerical breakdown";
      break;
    case AZ_loss:
      state = SolveCrashed;
      msg = "failed: numerical loss of precision";
      break;
    case AZ_ill_cond:
      state = SolveCrashed;
      msg = "failed: ill-conditioned Hessenberg matrix in GMRES";
      break;
    case AZ_maxits:
      state = SolveFailedToConverge;
      msg = "failed: maxiters reached without converged";
      break;
  }
  SolverState<double> rtn(state, "Aztec solver " + msg, (int) status[AZ_its],
    status[AZ_r]);
  return rtn;
}


void AztecSolver::setUserPrec(const LinearOperator<double>& P,
  const LinearSolver<double>& pSolver)
{
  if (useUserPrec_)
  {
    userPrec_ = rcp(new Epetra::Epetra_TSFOperator(P, pSolver));
  }
  else
  {
    TEST_FOR_EXCEPTION(!useUserPrec_, std::runtime_error,
      "Attempt to set user-defined preconditioner "
      "after another preconditioner has been specified");
  }
}

void AztecSolver::initParamMap()
{
  static bool first = true;
  if (first)
  {
    paramMap()["Method"]=AZ_solver;
    paramMap()["CG"]=AZ_cg;
    paramMap()["GMRES"]=AZ_gmres;
    paramMap()["CGS"]=AZ_cgs;
    paramMap()["TFQMR"]=AZ_tfqmr;
    paramMap()["BICGSTAB"]=AZ_bicgstab;
    paramMap()["Direct"]=AZ_lu;
    paramMap()["Precond"]=AZ_precond;
    paramMap()["None"]=AZ_none;
    paramMap()["Jacobi"]=AZ_Jacobi;
    paramMap()["Neumann Series"]=AZ_Neumann;
    paramMap()["Symmetric Gauss-Seidel"]=AZ_sym_GS;
    paramMap()["Least-Squares Polynomial"]=AZ_ls;
    paramMap()["Recursive Iterate"]=AZ_recursive_iterate;
    paramMap()["Domain Decomposition"]=AZ_dom_decomp;
    paramMap()["Subdomain Solver"]=AZ_subdomain_solve;
    paramMap()["Approximate Sparse LU"]=AZ_lu;
    paramMap()["Saad ILUT"]=AZ_ilut;
    paramMap()["ILU"]=AZ_ilu;
    paramMap()["RILU"]=AZ_rilu;
    paramMap()["Block ILU"]=AZ_bilu;
    paramMap()["Incomplete Cholesky"]=AZ_icc;
    paramMap()["Residual Scaling"]=AZ_conv;
    paramMap()["Initial"]=AZ_r0;
    paramMap()["RHS"]=AZ_rhs;
    paramMap()["Matrix"]=AZ_Anorm;
    paramMap()["Solution"]=AZ_sol;
    paramMap()["No Scaling"]=AZ_noscaled;
    paramMap()["Verbosity"]=AZ_output;
    paramMap()["All"]=AZ_all;
    paramMap()["Silent"]=AZ_none;
    paramMap()["Warnings"]=AZ_warnings;
    paramMap()["Final Residual"]=AZ_last;
    paramMap()["Graph Fill"]=AZ_graph_fill;
    paramMap()["Max Iterations"]=AZ_max_iter;
    paramMap()["Polynomial Order"]=AZ_poly_ord;
    paramMap()["Overlap"]=AZ_overlap;
    paramMap()["Overlap Type"]=AZ_type_overlap;
    paramMap()["Standard"]=AZ_standard;
    paramMap()["Symmetric"]=AZ_symmetric;
    paramMap()["Restart Size"]=AZ_kspace;
    paramMap()["Reorder ILU"]=AZ_reorder;
    paramMap()["Keep Factorization"]=AZ_keep_info;
    paramMap()["GMRES Orthogonalization"]=AZ_orthog;
    paramMap()["Classical Gram-Schmidt"]=AZ_classic;
    paramMap()["Modified Gram-Schmidt"]=AZ_modified;
    paramMap()["Auxiliary Vector"]=AZ_aux_vec;
    paramMap()["Residual"]=AZ_resid;
    paramMap()["Random"]=AZ_rand;
    paramMap()["Tolerance"]=AZ_tol;
    paramMap()["Drop Tolerance"]=AZ_drop;
    paramMap()["Fill Ratio"]=AZ_ilut_fill;
    paramMap()["Damping"]=AZ_omega;

    first = false;
  }
}


