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
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_mmio.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"

using namespace EpetraExt;
namespace EpetraExt {

int BlockMapToMatrixMarketFile( const char *filename, const Epetra_BlockMap & map, 
				 const char * mapName,
				 const char *mapDescription, 
				 bool writeHeader) {
  int M = map.NumGlobalElements();
  int N = 1;
  if (map.MaxElementSize()>1) N = 2; // Non-trivial block map, store element sizes in second column

  FILE * handle = 0;

  if (map.Comm().MyPID()==0) { // Only PE 0 does this section

    handle = fopen(filename,"w");
    if (!handle) return(-1);
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_integer(&matcode);

    if (writeHeader==true) { // Only write header if requested (true by default)
    
      if (mm_write_banner(handle, matcode)) return(-1);
      
      if (mapName!=0) fprintf(handle, "%% \n%% %s\n", mapName);
      if (mapDescription!=0) fprintf(handle, "%% %s\n%% \n", mapDescription);

    }
  }
    
  if (writeHeader==true) { // Only write header if requested (true by default)

    // Make an Epetra_IntVector of length numProc such that all elements are on PE 0 and
    // the ith element is NumMyElements from the ith PE

    Epetra_Map map1(-1, 1, 0, map.Comm()); // map with one element on each processor
    int length = 0;
    if (map.Comm().MyPID()==0) length = map.Comm().NumProc();
    Epetra_Map map2(-1, length, 0, map.Comm());
    Epetra_Import lengthImporter(map2, map1);
    Epetra_IntVector v1(map1);
    Epetra_IntVector v2(map2);
    v1[0] = map.NumMyElements();
    if (v2.Import(v1, lengthImporter, Insert)) return(-1);
    if (map.Comm().MyPID()==0) { 
      fprintf(handle, "%s", "%Format Version:\n");
      //int version = 1; // We may change the format scheme at a later date.
      fprintf(handle, "%% %d \n", map.Comm().NumProc());
      fprintf(handle, "%s", "%NumProc: Number of processors:\n");
      fprintf(handle, "%% %d \n", map.Comm().NumProc());
      fprintf(handle, "%s", "%MaxElementSize: Maximum element size:\n");
      fprintf(handle, "%% %d \n", map.MaxElementSize());
      fprintf(handle, "%s", "%MinElementSize: Minimum element size:\n");
      fprintf(handle, "%% %d \n", map.MinElementSize());
      fprintf(handle, "%s", "%IndexBase: Index base of map:\n");
      fprintf(handle, "%% %d \n", map.IndexBase());
      fprintf(handle, "%s", "%NumGlobalElements: Total number of GIDs in map:\n");
      fprintf(handle, "%% %d \n", map.NumGlobalElements());
      fprintf(handle, "%s", "%NumMyElements: BlockMap lengths per processor:\n");
      for ( int i=0; i< v2.MyLength(); i++) fprintf(handle, "%% %d\n", v2[i]);
      
      if (mm_write_mtx_array_size(handle, M, N)) return(-1);
    }
  }
  if (BlockMapToHandle(handle, map)) return(-1); // Everybody calls this routine
  
  if (map.Comm().MyPID()==0) // Only PE 0 opened a file
    if (fclose(handle)) return(-1);
  return(0);
}

int BlockMapToHandle(FILE * handle, const Epetra_BlockMap & map) {

  const Epetra_Comm & comm = map.Comm();
  int numProc = comm.NumProc();
  bool doSizes = !map.ConstantElementSize();

  if (numProc==1) {
    int * myElements = map.MyGlobalElements();
    int * elementSizeList = 0;
    if (doSizes) elementSizeList = map.ElementSizeList();
    return(writeBlockMap(handle, map.NumGlobalElements(), myElements, elementSizeList, doSizes));
  }

  int numRows = map.NumMyElements();
  
  Epetra_Map allGidsMap(-1, numRows, 0,comm);
  
  Epetra_IntVector allGids(allGidsMap);
  for (int i=0; i<numRows; i++) allGids[i] = map.GID(i);
  
  Epetra_IntVector allSizes(allGidsMap);
  for (int i=0; i<numRows; i++) allSizes[i] = map.ElementSize(i);
  
  // Now construct a Map on PE 0 by strip-mining the rows of the input matrix map.
  int numChunks = numProc;
  int stripSize = allGids.GlobalLength()/numChunks;
  int remainder = allGids.GlobalLength()%numChunks;
  int curStart = 0;
  int curStripSize = 0;
  Epetra_IntSerialDenseVector importGidList;
  Epetra_IntSerialDenseVector importSizeList;
  if (comm.MyPID()==0) {
    importGidList.Size(stripSize+1); // Set size of vector to max needed
    if (doSizes) importSizeList.Size(stripSize+1); // Set size of vector to max needed
  }
  for (int i=0; i<numChunks; i++) {
    if (comm.MyPID()==0) { // Only PE 0 does this part
      curStripSize = stripSize;
      if (i<remainder) curStripSize++; // handle leftovers
      for (int j=0; j<curStripSize; j++) importGidList[j] = j + curStart;
      curStart += curStripSize;
    }
    // The following import map will be non-trivial only on PE 0.
    Epetra_Map importGidMap(-1, curStripSize, importGidList.Values(), 0, comm);
    Epetra_Import gidImporter(importGidMap, allGidsMap);
    
    Epetra_IntVector importGids(importGidMap);
    if (importGids.Import(allGids, gidImporter, Insert)) return(-1); 
    Epetra_IntVector importSizes(importGidMap);
    if (doSizes) if (importSizes.Import(allSizes, gidImporter, Insert)) return(-1); 
    
    // importGids (and importSizes, if non-trivial block map)
    // now have a list of GIDs (and sizes, respectively) for the current strip of map.
    
    int * myElements = importGids.Values();
    int * elementSizeList = 0;
    if (doSizes) elementSizeList = importSizes.Values();
    // Finally we are ready to write this strip of the map to file
    writeBlockMap(handle, importGids.MyLength(), myElements, elementSizeList, doSizes);
  }
  return(0);
}
int writeBlockMap(FILE * handle, int length, const int * v1, const int * v2, bool doSizes) {

  for (int i=0; i<length; i++) {
    fprintf(handle, "%d", v1[i]);
    if (doSizes) fprintf(handle, " %d", v2[i]);
    fprintf(handle, "%s", "\n");
  }
  return(0);
}
} // namespace EpetraExt
