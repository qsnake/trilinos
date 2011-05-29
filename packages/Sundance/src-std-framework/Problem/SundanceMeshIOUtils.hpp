#ifndef SUNDANCE_MESHIO_UTILS_HPP
#define SUNDANCE_MESHIO_UTILS_HPP

#include "Sundance.hpp"
#include "SundanceFileIOChacoPartitioner.hpp"
#include "SundanceCellLIDMappedFieldWrapper.hpp"

namespace Sundance
{

/** Read the node-based fields from a mesher */
Expr readNodalFields(const MeshSource& mesher, const Mesh& mesh,
  const VectorType<double>& vecType);


/** Read a 2D field stored in simple ASCII format
 * 
 * nx ny                 # num x pts, num y pts
 * u1_1 u2_1 ... uN_1    # data at node 1
 * u1_2 u2_2 ... uN_2    # data at node 2
 * 
*/
Expr readSerialGridField(const std::string& gridFile, 
  double ax, double bx,
  double ay, double by,
  const VectorType<double>& vecType,
  const MeshType& meshType,
  Mesh& mesh);


/** This function reads in an exodus file, discretizes a function on 
 * it, and computes some functional on it. It then writes the field, 
 * reads it back, and computes the same functional on the new copy. 
 * The return value is the difference between the two functionals
 * which should be zero if all is running correctly. */
double readbackTester(const std::string& infile, const MPIComm& comm) ;


/** Partition a mesh and write the parts to exodus files */
void serialPartition(
  const RCP<SerialPartitionerBase>& part,
  int numProc,
  const MeshSource& mesher, 
  const std::string& outfile);
}

#endif
