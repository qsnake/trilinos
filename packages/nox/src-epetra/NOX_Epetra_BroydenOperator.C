//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "NOX_Epetra_BroydenOperator.H"

// EpetraExt includes for dumping Epetra objects
#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#endif
#endif

#include "Epetra_Map.h"

using namespace NOX;
using namespace NOX::Epetra;

BroydenOperator::BroydenOperator(
       Teuchos::ParameterList & nlParams_                ,
       const Teuchos::RCP<NOX::Utils>& utils_    ,
       Epetra_Vector & solnVec                           ,
       const Teuchos::RCP<Epetra_CrsMatrix>& mat ,
       bool verbose_ ) :
  verbose        ( verbose_                                  ) ,
  crsMatrix      ( Teuchos::rcp( new Epetra_CrsMatrix(*mat)) ) ,
  nlParams       ( nlParams_                                 ) ,
  utils          ( utils_                                    ) ,
  prePostOperator( utils, nlParams.sublist("Solver Options") ) ,
  label          ( "NOX::Epetra::BroydenOperator"            ) ,
  isValidStep    ( false                                     ) ,
  isValidYield   ( false                                     ) ,
  isValidBroyden ( false                                     ) ,
  entriesRemoved ( mat->NumMyRows(), false                   )
{
  initialize( nlParams, solnVec );
}

//-----------------------------------------------------------------------------

BroydenOperator::BroydenOperator(const BroydenOperator & bOp) :
  verbose        ( bOp.verbose        ) ,
  stepVec        ( bOp.stepVec        ) ,
  yieldVec       ( bOp.yieldVec       ) ,
  workVec        ( bOp.workVec        ) ,
  oldX           ( bOp.oldX           ) ,
  oldF           ( bOp.oldF           ) ,
  crsMatrix      ( bOp.crsMatrix      ) ,
  jacIntPtr      ( bOp.jacIntPtr      ) ,
  jacMatrixPtr   ( bOp.jacMatrixPtr   ) ,
  precIntPtr     ( bOp.precIntPtr     ) ,
  precMatrixPtr  ( bOp.precMatrixPtr  ) ,
  nlParams       ( bOp.nlParams       ) ,
  utils          ( bOp.utils          ) ,
  prePostOperator( utils, nlParams.sublist("Solver Options") ),
  label          ( "NOX::Epetra::BroydenOperator" ),
  isValidStep    ( bOp.isValidStep    ) ,
  isValidYield   ( bOp.isValidYield   ) ,
  isValidBroyden ( bOp.isValidBroyden ) ,
  entriesRemoved ( bOp.entriesRemoved )
{ 
}

//-----------------------------------------------------------------------------

bool
BroydenOperator::initialize( Teuchos::ParameterList & nlParams, const Epetra_Vector & vec )
{
  stepVec  = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  yieldVec = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  workVec  = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  oldX     = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  oldF     = Teuchos::rcp( new NOX::Epetra::Vector(vec) );

  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  // RPP 9/20/2005: This is a very bad idea!  It breaks the rcp and
  // user's expectations.  For now we will have to create rcp without
  // ownership.  What happens if a user write their own PPO?
  Teuchos::RCP<NOX::Abstract::PrePostOperator> me = Teuchos::rcp(this, false);
  nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator", me);

  return true;
}

//-----------------------------------------------------------------------------

BroydenOperator::~BroydenOperator()
{
}

const char* BroydenOperator::Label () const
{
  return label.c_str();
}

int BroydenOperator::SetUseTranspose(bool UseTranspose)
{
  return crsMatrix->SetUseTranspose(UseTranspose);
}

int BroydenOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return crsMatrix->Apply(X, Y);
}

int BroydenOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return crsMatrix->ApplyInverse(X, Y);
}

bool BroydenOperator::UseTranspose() const
{
  return crsMatrix->UseTranspose();
}

bool BroydenOperator::HasNormInf() const
{
  return crsMatrix->HasNormInf();
}

const Epetra_Map& BroydenOperator::OperatorDomainMap() const
{
  return crsMatrix->OperatorDomainMap();
}

const Epetra_Map& BroydenOperator::OperatorRangeMap() const
{
  return crsMatrix->OperatorRangeMap();
}

bool BroydenOperator::Filled() const
{
  return crsMatrix->Filled();
}

int BroydenOperator::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  return crsMatrix->NumMyRowEntries(MyRow, NumEntries);
}

int BroydenOperator::MaxNumEntries() const
{
  return crsMatrix->MaxNumEntries();
}

inline int BroydenOperator::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  return crsMatrix->ExtractMyRowCopy(MyRow, Length, NumEntries, Values, Indices);
}

int BroydenOperator::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  return crsMatrix->ExtractDiagonalCopy(Diagonal);
}

int BroydenOperator::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return crsMatrix->Multiply(TransA, X, Y);
}

int BroydenOperator::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,  Epetra_MultiVector& Y) const
{
  return crsMatrix->Solve(Upper, Trans, UnitDiagonal, X, Y);
}

int BroydenOperator::InvRowSums(Epetra_Vector& x) const
{
  return crsMatrix->InvRowSums(x);
}

int BroydenOperator::LeftScale(const Epetra_Vector& x)
{
  return crsMatrix->LeftScale(x);
}

int BroydenOperator::InvColSums(Epetra_Vector& x) const
{
  return crsMatrix->InvColSums(x);
}

int BroydenOperator::RightScale(const Epetra_Vector& x)
{
  return crsMatrix->RightScale(x);
}

double BroydenOperator::NormInf() const
{
  return crsMatrix->NormInf();
}

double BroydenOperator::NormOne() const
{
  return crsMatrix->NormOne();
}

int BroydenOperator::NumGlobalNonzeros() const
{
  return crsMatrix->NumGlobalNonzeros();
}

int BroydenOperator::NumGlobalRows() const
{
  return crsMatrix->NumGlobalRows();
}

int BroydenOperator::NumGlobalCols() const
{
  return crsMatrix->NumGlobalCols();
}

int BroydenOperator::NumGlobalDiagonals() const
{
  return crsMatrix->NumGlobalDiagonals();
}

int BroydenOperator::NumMyNonzeros() const
{
  return crsMatrix->NumMyNonzeros();
}

int BroydenOperator::NumMyRows() const
{
  return crsMatrix->NumMyRows();
}

int BroydenOperator::NumMyCols() const
{
  return crsMatrix->NumMyCols();
}

int BroydenOperator::NumMyDiagonals() const
{
  return crsMatrix->NumMyDiagonals();
}

bool BroydenOperator::LowerTriangular() const
{
  return crsMatrix->LowerTriangular();
}

bool BroydenOperator::UpperTriangular() const
{
  return crsMatrix->UpperTriangular();
}

const Epetra_Comm& BroydenOperator::Comm() const
{
  return crsMatrix->Comm();
}

const Epetra_Map& BroydenOperator::RowMatrixRowMap() const
{
  return crsMatrix->RowMatrixRowMap();
}

const Epetra_Map& BroydenOperator::RowMatrixColMap() const
{
  return crsMatrix->RowMatrixColMap();
}

const Epetra_Import* BroydenOperator::RowMatrixImporter() const
{
  return crsMatrix->RowMatrixImporter();
}

const Epetra_BlockMap& BroydenOperator::Map() const
{
  return crsMatrix->Map();
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computeJacobian( const Epetra_Vector & x, Epetra_Operator& Jac )
{
  bool ok = false;

  ok = computeSparseBroydenUpdate();

  vector<ReplacementInterface *>::iterator iter = replacementInterfaces.begin() ,
                                      iter_end  = replacementInterfaces.end()     ;

  for( ; iter_end != iter; ++iter )
  {
    Teuchos::RCP<const Epetra_CrsMatrix> pMat = (*iter)->getReplacementValuesMatrix(x, ReplacementInterface::JACOBIAN );
    replaceBroydenMatrixValues( *pMat );
  }

  return ok;

}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computePreconditioner( const Epetra_Vector & x,
					Epetra_Operator& Prec,
                                        Teuchos::ParameterList * pList )
{
  bool ok = false;

  ok = computeSparseBroydenUpdate();

  vector<ReplacementInterface *>::iterator iter = replacementInterfaces.begin() ,
                                      iter_end  = replacementInterfaces.end()     ;

  for( ; iter_end != iter; ++iter )
  {
    Teuchos::RCP<const Epetra_CrsMatrix> pMat = (*iter)->getReplacementValuesMatrix(x, ReplacementInterface::PRECONDITIONER );
    replaceBroydenMatrixValues( *pMat );
  }

  return ok;

}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPreSolve( const NOX::Solver::Generic & solver)
{
  // Reset valid step and yield flags. NOTE: we leave the Broyden valid flag
  // alone as it may be valid from a previous solve
  isValidStep  = false;
  isValidYield = false;

  prePostOperator.runPreSolve( solver );

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPreIterate( const NOX::Solver::Generic & solver)
{
  prePostOperator.runPreIterate( solver );

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPostIterate( const NOX::Solver::Generic & solver)
{
  // Get and update using the solver object.

  if( solver.getNumIterations() > 0 )
  {
    // To consistently compute changes to the step and the yield, eg effects of
    // linesearch or other globalizations, we need to get and subtract the old and
    // new solution states and corresponding residuals 

    const Abstract::Group & oldSolnGrp = solver.getPreviousSolutionGroup();
    const Abstract::Group & solnGrp    = solver.getSolutionGroup();

    // Set the Step vector
    *workVec = solnGrp.getX();
    workVec->update(-1.0, oldSolnGrp.getX(), 1.0);

    setStepVector( *workVec );

    // Set the Yield vector
    *workVec = solnGrp.getF();
    workVec->update(-1.0, oldSolnGrp.getF(), 1.0);

    setYieldVector( *workVec );

    // Use of these updated vectors occurs when computeJacobian or computePreconditioner 
    // gets called
  }

  prePostOperator.runPostIterate( solver );

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPostSolve( const NOX::Solver::Generic & solver)
{
  prePostOperator.runPostSolve( solver );

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setStepVector( Epetra_Vector & vec )
{
  stepVec->getEpetraVector() = vec;
  isValidStep    = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setStepVector( NOX::Epetra::Vector & vec )
{
  *stepVec = vec;
  isValidStep    = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setYieldVector( NOX::Epetra::Vector & vec )
{
  *yieldVec = vec;
  isValidYield   = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setYieldVector( Epetra_Vector & vec )
{
  yieldVec->getEpetraVector() = vec;
  isValidYield   = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computeSparseBroydenUpdate()
{

  if( isValidBroyden )
    return true;
  else if( isValidStep && isValidYield )
  {
    // Do the Broyden update to our matrix
    int ierr = crsMatrix->Multiply( false, stepVec->getEpetraVector(), workVec->getEpetraVector() );
    if( ierr )
    {
      cout << "ERROR: NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate(...) "
           << "- crsMatrix->Multiply() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }

    int      numEntries  = crsMatrix->NumMyCols();
    int    * indices ;
    double * values  ;

    for( int row = 0; row < crsMatrix->NumMyRows(); ++row )
    {
      ierr = crsMatrix->ExtractMyRowView( row, numEntries, values, indices );
      if( ierr )
      {
        cout << "ERROR (" << ierr << ") : NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate() "
             << "- crsMatrix->ExtractGlobalRowView(...) failed for row --> " << row << endl;
        throw "NOX Error: Broyden Update Failed";
      }

      double diffVal = yieldVec->getEpetraVector()[row] - workVec->getEpetraVector()[row];

      double rowStepInnerProduct = 0.0;

      if( entriesRemoved[row] )
      {
        list<int>::iterator     iter = retainedEntries[row].begin() ,
                            iter_end = retainedEntries[row].end()    ;

        for( ; iter_end != iter; ++iter )
          rowStepInnerProduct += 
            stepVec->getEpetraVector()[indices[*iter]] * stepVec->getEpetraVector()[indices[*iter]];

        for( iter = retainedEntries[row].begin(); iter_end != iter; ++iter )
          values[*iter] += diffVal * stepVec->getEpetraVector()[indices[*iter]] / rowStepInnerProduct;
      }
      else
      {
        for( int col = 0; col < numEntries; ++col )
          rowStepInnerProduct += 
            stepVec->getEpetraVector()[indices[col]] * stepVec->getEpetraVector()[indices[col]];

        for( int col = 0; col < numEntries; ++col )
          (*values++) += diffVal * stepVec->getEpetraVector()[(*indices++)] / rowStepInnerProduct;
      }
    }

    // Our %Broyden matrix has been updated and is now ready to use as a 
    // preconditioning matrix or as the Jacobian.

    if( verbose )
      cout << "\t... BroydenOperator::computeSparseBroydenUpdate()..." << endl;

  }
  else
  {
    cout << "Warning: NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate(...) "
         << "- either the step vector or the yield vector or both is not valid." << endl;
    cout <<  "Leaving existing matrix unchanged." << endl;
  }

  // Use EpetraExt to dump linear system if debuggging
#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT

  static int broydenOutputCount;

  Teuchos::ParameterList & broydenParamas = 
    nlParams.sublist("Direction").sublist("Newton").sublist("Broyden Op");

  if( broydenParamas.get("Write Broyden Info", false) ) 
  {
    std::ostringstream outputNumber;
    outputNumber << broydenOutputCount++ ;
    
    std::string prefixName = broydenParamas.get("Write Broyden Info File Prefix", "BroydenOp");
    std::string postfixName = outputNumber.str();
    postfixName += ".mm";

    std::string mapFileName = prefixName + "_Map_"      + postfixName;
    std::string matFileName = prefixName + "_Matrix_"   + postfixName;    
    std::string dxFileName  = prefixName + "_dX_"       + postfixName;
    std::string dfFileName  = prefixName + "_dF_"       + postfixName;

    Epetra_RowMatrix * printMatrix = NULL;
    printMatrix = dynamic_cast<Epetra_RowMatrix*>(crsMatrix.get()); 

    if( NULL == printMatrix ) 
    {
      cout << "Error: NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate() - "
	   << "Could not get a valid crsMatrix!\n"
	   << "Please set the \"Write Linear System\" parameter to false." << endl;
      throw "NOX Error";
    }

    EpetraExt::BlockMapToMatrixMarketFile(mapFileName.c_str(), 
					  printMatrix->RowMatrixRowMap()); 
    EpetraExt::RowMatrixToMatrixMarketFile(matFileName.c_str(), *printMatrix, 
					   "test matrix", "Broyden Matrix XXX");
    EpetraExt::MultiVectorToMatrixMarketFile(dxFileName.c_str(), 
					     stepVec->getEpetraVector());
    EpetraExt::MultiVectorToMatrixMarketFile(dfFileName.c_str(), 
					     yieldVec->getEpetraVector());

  }
#endif
#endif

  return true;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::resetBroydenMatrix( const Epetra_CrsMatrix & mat )
{
  *crsMatrix = mat;

  return;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::isStep() const
{
  return isValidStep;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::isYield() const
{
  return isValidYield;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::isBroyden() const
{
  return isValidBroyden;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::removeEntriesFromBroydenUpdate( const Epetra_CrsGraph & graph )
{

  int numRemoveIndices ;
  int * removeIndPtr   ;
  int ierr             ;

  cout << graph << endl;

  for( int row = 0; row < graph.NumMyRows(); ++row) 
  {
    ierr = graph.ExtractMyRowView( row, numRemoveIndices, removeIndPtr );
    if( ierr )
    {
      cout << "ERROR (" << ierr << ") : "
           << "NOX::Epetra::BroydenOperator::removeEntriesFromBroydenUpdate(...)"
           << " - Extract indices error for row --> " << row << endl;
      throw "NOX Broyden Operator Error";
    }

    if( 0 != numRemoveIndices )
    {
      // Create a map for quick queries
      map<int, bool> removeIndTable;
      for( int k = 0; k < numRemoveIndices; ++k )
        removeIndTable[ graph.ColMap().GID(removeIndPtr[k]) ] = true;

      // Get our matrix column indices for the current row
      int numOrigIndices = 0;
      int * indPtr;

      ierr = crsMatrix->Graph().ExtractMyRowView( row, numOrigIndices, indPtr );
      if( ierr )
      {
        cout << "ERROR (" << ierr << ") : "
             << "NOX::Epetra::BroydenOperator::removeEntriesFromBroydenUpdate(...)"
             << " - Extract indices error for row --> " << row << endl;
        throw "NOX Broyden Operator Error";
      }

      // Remove appropriate active entities
      if( retainedEntries.end() == retainedEntries.find(row) )
      {
        list<int> inds;

        for( int k = 0; k < numOrigIndices; ++k )
        {
          if( removeIndTable.end() == removeIndTable.find( crsMatrix->Graph().ColMap().GID(indPtr[k]) ) )
            inds.push_back(k);
        }

        retainedEntries[row] = inds;
      }
      else
      {
        list<int> & inds = retainedEntries[row];

        list<int>::iterator iter     = inds.begin() ,
                            iter_end = inds.end()    ;

        for( ; iter_end != iter; ++iter )
        {
          if( !removeIndTable[ *iter ] )
            inds.remove( *iter );
        }
      }

      entriesRemoved[row] = true;
    }
  }

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::replaceBroydenMatrixValues( const Epetra_CrsMatrix & mat)
{
  double * values    ;
  int    * indices   ;
  int     numEntries ;
  int     ierr       ;

  for( int row = 0; row < mat.NumMyRows(); ++row) 
  {
    ierr = mat.ExtractMyRowView(row, numEntries, values, indices);
    ierr += crsMatrix->ReplaceGlobalValues(row, numEntries, values, indices);
    if( ierr )
    {
      cout << "ERROR (" << ierr << ") : "
           << "NOX::Epetra::BroydenOperator::replaceBroydenMatrixValues(...)"
           << " - Extract or Replace values error for row --> "
           << row << endl;
      throw "NOX Broyden Operator Error";
    }
  }
}

//-----------------------------------------------------------------------------

#ifdef HAVE_NOX_DEBUG

void 
BroydenOperator::outputActiveEntries()
{

  const Epetra_Comm & comm = stepVec->getEpetraVector().Comm();
  int numProcs = comm.NumProc();
  int myPID = nlParams.sublist("Printing").get("MyPID", 0);

  map<int, list<int> >::iterator iter     = retainedEntries.begin(),
                                 iter_end = retainedEntries.end()    ;

  for( int pid = 0; pid < numProcs; ++pid )
  {
    if( pid == myPID )
    {
      for( ; iter_end != iter; ++iter )
      {
        int row = (*iter).first;

        // Get our matrix column indices for the current row
        int numIndices ;
        int * indPtr   ;

        crsMatrix->Graph().ExtractMyRowView( row, numIndices, indPtr );

        list<int> & entries = (*iter).second;

        list<int>::iterator colIter     = entries.begin() ,
                            colIter_end = entries.end()    ;

        for( ; colIter_end != colIter; ++colIter )
          cout << "[" << pid << "] \t " << row << " \t " << *colIter << " ("
               << indPtr[*colIter] << ")" << endl;
      }
    }
    comm.Barrier();
  }

  return;
}

//-----------------------------------------------------------------------------

#endif
