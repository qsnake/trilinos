//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
//@HEADER
#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_mmio.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"

using namespace EpetraExt;
namespace EpetraExt {

int OperatorToMatlabFile( const char *filename, const Epetra_Operator & A) {

  // Simple wrapper to make it clear what can be used to write to Matlab format
  EPETRA_CHK_ERR(OperatorToMatrixMarketFile(filename, A, 0, 0, false));
  return(0);
}

int OperatorToMatrixMarketFile( const char *filename, const Epetra_Operator & A, 
				 const char * matrixName,
				 const char *matrixDescription, 
				 bool writeHeader) {

  const Epetra_Map & domainMap = A.OperatorDomainMap();
  const Epetra_Map & rangeMap = A.OperatorRangeMap();

  if (!domainMap.UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
  if (!rangeMap.UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
  
  int M = rangeMap.NumGlobalElements();
  int N = domainMap.NumGlobalElements();

  FILE * handle = 0;

  // To get count of nonzero terms we do multiplies ...
  int nz = 0;
  if (get_nz(A, nz)) {EPETRA_CHK_ERR(-1);}

  if (domainMap.Comm().MyPID()==0) { // Only PE 0 does this section

    handle = fopen(filename,"w");
    if (!handle) {EPETRA_CHK_ERR(-1);}
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    

    if (writeHeader==true) { // Only write header if requested (true by default)
    
      if (mm_write_banner(handle, matcode)!=0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
      
      if (matrixName!=0) fprintf(handle, "%% \n%% %s\n", matrixName);
      if (matrixDescription!=0) fprintf(handle, "%% %s\n%% \n", matrixDescription);
      
      if (mm_write_mtx_crd_size(handle, M, N, nz)!=0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
    }
  }
    
  if (OperatorToHandle(handle, A)!=0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}// Everybody calls this routine

  if (handle!=0) fclose(handle);
  return(0);
}

int OperatorToHandle(FILE * handle, const Epetra_Operator & A) {

  const Epetra_Map & domainMap = A.OperatorDomainMap();
  const Epetra_Map & rangeMap = A.OperatorRangeMap();
  int N = domainMap.NumGlobalElements();

  //cout << "rangeMap = " << rangeMap << endl;
  Epetra_Map rootRangeMap = Epetra_Util::Create_Root_Map(rangeMap);
  //cout << "rootRangeMap = " << rootRangeMap << endl;
  Epetra_Map rootDomainMap = Epetra_Util::Create_Root_Map(domainMap, -1); // Replicate on all processors
  Epetra_Import importer(rootRangeMap, rangeMap);

  int chunksize = 5; // Let's do multiple RHS at a time
  int numchunks = N/chunksize;
  int rem = N%chunksize;

  if (rem>0) {
    Epetra_MultiVector xrem(domainMap, rem);
    Epetra_MultiVector yrem(rangeMap, rem);
    Epetra_MultiVector yrem1(rootRangeMap, rem);
    // Put 1's in slots
    for (int j=0; j<rem; j++) {
      int curGlobalCol = rootDomainMap.GID(j); // Should return same value on all processors
      if (domainMap.MyGID(curGlobalCol)) {
	int curCol = domainMap.LID(curGlobalCol);
	xrem[j][curCol] = 1.0;
      }
    }
    EPETRA_CHK_ERR(A.Apply(xrem, yrem));
    EPETRA_CHK_ERR(yrem1.Import(yrem, importer, Insert));
    EPETRA_CHK_ERR(writeOperatorStrip(handle, yrem1, rootDomainMap, rootRangeMap, 0));
  }

  if (numchunks>0) {
    Epetra_MultiVector x(domainMap, chunksize);
    Epetra_MultiVector y(rangeMap, chunksize);
    Epetra_MultiVector y1(rootRangeMap, chunksize);
    for (int ichunk = 0; ichunk<numchunks; ichunk++) {
      int startCol = ichunk*chunksize+rem;
      // Put 1's in slots
      for (int j=0; j<chunksize; j++) {
	int curGlobalCol = rootDomainMap.GID(startCol+j); // Should return same value on all processors
	if (domainMap.MyGID(curGlobalCol)){
	  int curCol = domainMap.LID(curGlobalCol);
	  x[j][curCol] = 1.0;
	}
      }
      EPETRA_CHK_ERR(A.Apply(x, y));
      EPETRA_CHK_ERR(y1.Import(y, importer, Insert));
      EPETRA_CHK_ERR(writeOperatorStrip(handle, y1, rootDomainMap, rootRangeMap, startCol));
      // Put 0's in slots
      for (int j=0; j<chunksize; j++) {
	int curGlobalCol = rootDomainMap.GID(startCol+j); // Should return same value on all processors
	if (domainMap.MyGID(curGlobalCol)){
	  int curCol = domainMap.LID(curGlobalCol);
	  x[j][curCol] = 0.0;
	}
      }
    }
  }

  return(0);
}
int writeOperatorStrip(FILE * handle, const Epetra_MultiVector & y, const Epetra_Map & rootDomainMap, const Epetra_Map & rootRangeMap, int startColumn) {

  int numRows = y.GlobalLength();
  int numCols = y.NumVectors();
  int ioffset = 1 - rootRangeMap.IndexBase(); // Matlab indices start at 1
  int joffset = 1 - rootDomainMap.IndexBase(); // Matlab indices start at 1
  if (y.Comm().MyPID()!=0) {
    if (y.MyLength()!=0) {EPETRA_CHK_ERR(-1);}
  }
  else {
    if (numRows!=y.MyLength()) {EPETRA_CHK_ERR(-1);}
    for (int j=0; j<numCols; j++) {
      int J = rootDomainMap.GID(j + startColumn) + joffset;
      for (int i=0; i<numRows; i++) {
	double val = y[j][i];
	if (val!=0.0) {
	  int I = rootRangeMap.GID(i) + ioffset;
	  fprintf(handle, "%d %d %22.16e\n", I, J, val);
	}
      }
    }
  }
  return(0);
}
int get_nz(const Epetra_Operator & A, int & nz) {
  
  const Epetra_Map & domainMap = A.OperatorDomainMap();
  const Epetra_Map & rangeMap = A.OperatorRangeMap();
    

  int N = domainMap.NumGlobalElements();
  Epetra_Map rootDomainMap = Epetra_Util::Create_Root_Map(domainMap, -1); // Replicate on all processors


  int chunksize = 5; // Let's do multiple RHS at a time
  int numchunks = N/chunksize;
  int rem = N%chunksize;

  int lnz = 0;
  if (rem>0) {
    Epetra_MultiVector xrem(domainMap, rem);
    Epetra_MultiVector yrem(rangeMap, rem);
    // Put 1's in slots
    for (int j=0; j<rem; j++) {
      int curGlobalCol = rootDomainMap.GID(j);
      if (domainMap.MyGID(curGlobalCol)) xrem[j][domainMap.LID(curGlobalCol)] = 1.0;
    }
    EPETRA_CHK_ERR(A.Apply(xrem, yrem));
    for (int j=0; j<rem; j++) {
      int mylength = yrem.MyLength();
      for (int i=0; i<mylength; i++) 
	if (yrem[j][i]!=0.0) lnz++;
    }
  }

  if (numchunks>0) {
    Epetra_MultiVector x(domainMap, chunksize);
    Epetra_MultiVector y(rangeMap, chunksize);
    for (int ichunk = 0; ichunk<numchunks; ichunk++) {
      int startCol = ichunk*chunksize+rem;
      // Put 1's in slots
      for (int j=0; j<chunksize; j++) {
	int curGlobalCol = rootDomainMap.GID(startCol+j);
	if (domainMap.MyGID(curGlobalCol)) x[j][domainMap.LID(curGlobalCol)] = 1.0;
      }
      EPETRA_CHK_ERR(A.Apply(x, y));
      for (int j=0; j<chunksize; j++) {
	int mylength = y.MyLength();
	for (int i=0; i<mylength; i++) 
	  if (y[j][i]!=0.0) lnz++;
      }
      // Put 0's in slots
      for (int j=0; j<chunksize; j++) {
	int curGlobalCol = rootDomainMap.GID(startCol+j);
	if (domainMap.MyGID(curGlobalCol)) x[j][domainMap.LID(curGlobalCol)] = 0.0;
      }
    }
  }
    
  // Sum up nonzero counts
  EPETRA_CHK_ERR(A.Comm().SumAll(&lnz, &nz, 1));
  return(0);
}
} // namespace EpetraExt
