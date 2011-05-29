
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#
# NONE


#
# B) Define the header and source files (and directories)
#

#
# src
#


SET(HEADERS "")
SET(SOURCES "")


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})



#
# Vectors headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-utils/Utilities)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Vectors)


APPEND_SET(HEADERS)

APPEND_SET(SOURCES)



APPEND_SET(HEADERS
  Vectors/RTOpPack_SUNDIALS_Helpers.hpp
  Vectors/RTOpPack_SUNDIALS_Ops.hpp
  Vectors/Thyra_SUNDIALS_Ops.hpp
  Vectors/TSFBCPartitionedVSBuilder.hpp
  Vectors/TSFGhostImporter.hpp
  Vectors/TSFGhostView.hpp
  Vectors/TSFIndexableVector.hpp
  Vectors/TSFLoadableBlockVector.hpp
  Vectors/TSFLoadableVector.hpp
  Vectors/TSFLinearCombinationDecl.hpp
  Vectors/TSFLinearCombinationImpl.hpp
  Vectors/TSFAccessibleVector.hpp
  Vectors/TSFProductVectorSpaceDecl.hpp
  Vectors/TSFProductVectorSpaceImpl.hpp
  Vectors/TSFRawDataAccessibleVector.hpp
  Vectors/TSFSequentialIteratorDecl.hpp
  Vectors/TSFSequentialIteratorImpl.hpp
  Vectors/TSFVectorDecl.hpp
  Vectors/TSFVectorImpl.hpp
  Vectors/TSFVectorOpsDecl.hpp
  Vectors/TSFVectorOpsImpl.hpp
  Vectors/TSFVectorSpaceDecl.hpp
  Vectors/TSFVectorSpaceImpl.hpp
  Vectors/TSFVectorType.hpp
  Vectors/TSFVectorTypeExtensions.hpp
  )


APPEND_SET(SOURCES
  Vectors/TSFLinearCombination.cpp
  Vectors/TSFProductVectorSpace.cpp
  Vectors/TSFSequentialIterator.cpp
  Vectors/TSFVector.cpp
  Vectors/TSFVectorOps.cpp
  Vectors/TSFVectorSpace.cpp
)



#
# Concrete files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Concrete)

APPEND_SET(HEADERS
  Concrete/EpetraTSFOperator.hpp
  Concrete/TSFDenseSerialMatrix.hpp
  Concrete/TSFDenseSerialMatrixFactory.hpp
  Concrete/TSFDenseSerialVector.hpp
  Concrete/TSFEpetraGhostImporter.hpp
  Concrete/TSFEpetraGhostView.hpp
  Concrete/TSFEpetraMatrixFactory.hpp
  Concrete/TSFEpetraMatrix.hpp
  Concrete/TSFEpetraMatrixMatrixProduct.hpp
  Concrete/TSFEpetraMatrixMatrixSum.hpp
  Concrete/TSFEpetraMatrixOps.hpp
  Concrete/TSFEpetraVector.hpp
  Concrete/TSFEpetraVectorSpace.hpp
  Concrete/TSFEpetraVectorType.hpp
  Concrete/TSFIfpackOperator.hpp
  Concrete/TSFMLOperator.hpp
  Concrete/TSFLAPACKGeneralMatrix.hpp
  Concrete/TSFSerialVector.hpp
  Concrete/TSFSerialVectorSpace.hpp
  Concrete/TSFSerialGhostView.hpp
  Concrete/TSFSerialGhostImporter.hpp
  Concrete/TSFSerialVectorType.hpp
  Concrete/TSFVectorSpace2EpetraMap.hpp
  )

APPEND_SET(SOURCES
  Concrete/EpetraTSFOperator.cpp
  Concrete/TSFDenseSerialVector.cpp
  Concrete/TSFDenseSerialMatrix.cpp
  Concrete/TSFDenseSerialMatrixFactory.cpp
  Concrete/TSFEpetraGhostImporter.cpp
  Concrete/TSFEpetraGhostView.cpp
  Concrete/TSFEpetraMatrix.cpp
  Concrete/TSFEpetraMatrixFactory.cpp
  Concrete/TSFEpetraMatrixMatrixProduct.cpp
  Concrete/TSFEpetraMatrixMatrixSum.cpp
  Concrete/TSFEpetraMatrixOps.cpp
  Concrete/TSFEpetraVector.cpp
  Concrete/TSFEpetraVectorSpace.cpp
  Concrete/TSFEpetraVectorType.cpp
  Concrete/TSFIfpackOperator.cpp
  Concrete/TSFMLOperator.cpp
  Concrete/TSFLAPACKGeneralMatrix.cpp
  Concrete/TSFSerialVector.cpp
  Concrete/TSFSerialVectorSpace.cpp
  Concrete/TSFSerialGhostImporter.cpp
  Concrete/TSFSerialVectorType.cpp
  Concrete/TSFVectorSpace2EpetraMap.cpp
  )

#
# Operators files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Operators)

APPEND_SET(HEADERS
  Operators/TSFBlockOperatorBaseDecl.hpp
  Operators/TSFCollectivelyConfigurableMatrixFactory.hpp
  Operators/TSFEuclideanOpWithBackwardsCompatibleApply.hpp
  Operators/TSFExplicitlyTransposeableOp.hpp
  Operators/TSFExplicitlyTransposedLTIProblemFactory.hpp
  Operators/TSFHomogeneouslyBlockedLinearOp.hpp
  Operators/TSFILUFactorizableOp.hpp
  Operators/TSFIncrementallyConfigurableMatrixFactory.hpp
  Operators/TSFInverseOperatorDecl.hpp
  Operators/TSFInverseOperatorImpl.hpp
  Operators/TSFInverseLTIOp.hpp
  Operators/TSFLinearOperatorDecl.hpp
  Operators/TSFLinearOperatorImpl.hpp
  Operators/TSFLoadableBlockOperatorDecl.hpp
  Operators/TSFLoadableBlockOperatorImpl.hpp
  Operators/TSFLoadableMatrix.hpp
  Operators/TSFLTIProblemFactoryBase.hpp
  Operators/TSFDefaultLTIProblemFactory.hpp
  Operators/TSFMatrixFactory.hpp
  Operators/TSFMultiVectorOperatorDecl.hpp
  Operators/TSFMultiVectorOperatorImpl.hpp
  Operators/TSFNonlinearOperator.hpp
  Operators/TSFNonlinearOperatorBase.hpp
  Operators/TSFNonmemberOpHelpersDecl.hpp
  Operators/TSFNonmemberOpHelpersImpl.hpp
  Operators/TSFOpWithBackwardsCompatibleApply.hpp
  Operators/TSFPartitionedMatrixFactory.hpp
  Operators/TSFPartitionedToMonolithicConverter.hpp
  Operators/TSFRowAccessibleOp.hpp
  Operators/TSFSimpleAddedOpDecl.hpp
  Operators/TSFSimpleAddedOpImpl.hpp
  Operators/TSFSimpleBlockOpDecl.hpp
  Operators/TSFSimpleBlockOpImpl.hpp
  Operators/TSFSimpleComposedOpDecl.hpp
  Operators/TSFSimpleComposedOpImpl.hpp
  Operators/TSFSimpleDiagonalOpDecl.hpp
  Operators/TSFSimpleDiagonalOpImpl.hpp
  Operators/TSFSimpleIdentityOpDecl.hpp
  Operators/TSFSimpleIdentityOpImpl.hpp
  Operators/TSFSimpleScaledOpDecl.hpp
  Operators/TSFSimpleScaledOpImpl.hpp
  Operators/TSFSimpleTransposedOpDecl.hpp
  Operators/TSFSimpleTransposedOpImpl.hpp
  Operators/TSFSimpleZeroOpDecl.hpp
  Operators/TSFSimpleZeroOpImpl.hpp
  Operators/TSFSimplifiedLinearOpBaseDecl.hpp 
  Operators/TSFSimplifiedLinearOpBaseImpl.hpp 
  )

APPEND_SET(SOURCES
  Operators/TSFInverseOperator.cpp
  Operators/TSFLinearOperator.cpp
  Operators/TSFLoadableBlockOperator.cpp
  Operators/TSFMultiVectorOperator.cpp
  Operators/TSFPartitionedMatrixFactory.cpp
  Operators/TSFSimpleAddedOp.cpp
  Operators/TSFSimpleBlockOp.cpp
  Operators/TSFSimpleComposedOp.cpp
  Operators/TSFSimpleDiagonalOp.cpp
  Operators/TSFSimpleIdentityOp.cpp
  Operators/TSFSimpleScaledOp.cpp
  Operators/TSFSimpleTransposedOp.cpp
  Operators/TSFSimpleZeroOp.cpp
  Operators/TSFSimplifiedLinearOpBase.cpp 
  )

#
# LinearSolvers files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/LinearSolvers)

APPEND_SET(HEADERS
  LinearSolvers/TSFAmesosSolver.hpp
  LinearSolvers/TSFAnasaziEigensolverDecl.hpp
  LinearSolvers/TSFAnasaziEigensolverImpl.hpp
  LinearSolvers/TSFAztecSolver.hpp
  LinearSolvers/TSFAztecSolver.hpp
  LinearSolvers/TSFBelosSolver.hpp
  LinearSolvers/TSFBICGSTABSolverDecl.hpp
  LinearSolvers/TSFBICGSTABSolverImpl.hpp
  LinearSolvers/TSFBlockTriangularSolverDecl.hpp
  LinearSolvers/TSFBlockTriangularSolverImpl.hpp
  LinearSolvers/TSFDenseLUSolver.hpp
  LinearSolvers/TSFEigensolver.hpp
  LinearSolvers/TSFEigensolverBase.hpp
  LinearSolvers/TSFIterativeSolver.hpp
  LinearSolvers/TSFKrylovSolver.hpp
  LinearSolvers/TSFLinearSolverDecl.hpp
  LinearSolvers/TSFLinearSolverImpl.hpp
  LinearSolvers/TSFLinearSolverBaseDecl.hpp
  LinearSolvers/TSFLinearSolverBaseImpl.hpp
  LinearSolvers/TSFLinearSolverBuilder.hpp
  LinearSolvers/TSFSolverState.hpp 
  )

APPEND_SET(SOURCES
  LinearSolvers/TSFAmesosSolver.cpp
  LinearSolvers/TSFAnasaziEigensolver.cpp
  LinearSolvers/TSFAztecSolver.cpp
  LinearSolvers/TSFBelosSolver.cpp
  LinearSolvers/TSFBICGSTABSolver.cpp
  LinearSolvers/TSFBlockTriangularSolver.cpp
  LinearSolvers/TSFDenseLUSolver.cpp
  LinearSolvers/TSFLinearSolver.cpp
  LinearSolvers/TSFLinearSolverBase.cpp
  LinearSolvers/TSFLinearSolverBuilder.cpp
  )

#
# Non-LinearSolvers files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/NonlinearSolvers)

APPEND_SET(HEADERS
  NonlinearSolvers/NOX_StatusTest_SafeCombo.H
  NonlinearSolvers/NOX_TSF_Group.H
  NonlinearSolvers/NOX_TSF.H
  NonlinearSolvers/NOX_TSF_StatusTestBuilder.H
  NonlinearSolvers/NOX_TSF_Vector.H
  NonlinearSolvers/TSFNOXSolver.H 
  )

APPEND_SET(SOURCES
  NonlinearSolvers/NOX_StatusTest_SafeCombo.C
  NonlinearSolvers/NOX_TSF_Group.C
  NonlinearSolvers/NOX_TSF_StatusTestBuilder.C
  NonlinearSolvers/NOX_TSF_Vector.C
  NonlinearSolvers/TSFNOXSolver.C       
  )

#
# testing files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Testing)

APPEND_SET(HEADERS
  Testing/Thyra_TestSpecifier.hpp
  Testing/TSFCompoundTester.hpp
  Testing/TSFGlobalAnd.hpp
  Testing/TSFHeatOperator1D.hpp
  Testing/TSFLinearCombinationTester.hpp
  Testing/TSFMatrixLaplacian1D.hpp
  Testing/TSFMatrixMatrixTester.hpp
  Testing/TSFOperatorBuilder.hpp
  Testing/TSFPoissonBoltzmannJacobian.hpp
  Testing/TSFPoissonBoltzmannOp.hpp
  Testing/TSFRandomBlockMatrixBuilderDecl.hpp
  Testing/TSFRandomBlockMatrixBuilderImpl.hpp
  Testing/TSFRandomSparseMatrixBuilderDecl.hpp
  Testing/TSFRandomSparseMatrixBuilderImpl.hpp
  Testing/TSFTesterBase.hpp
  Testing/TSFVectorTester.hpp 
  )

APPEND_SET(SOURCES
  Testing/TSFGlobalAnd.cpp
  Testing/TSFHeatOperator1D.cpp
  Testing/TSFMatrixLaplacian1D.cpp
  Testing/TSFPoissonBoltzmannJacobian.cpp
  Testing/TSFPoissonBoltzmannOp.cpp 
  Testing/TSFRandomBlockMatrixBuilder.cpp
  Testing/TSFRandomSparseMatrixBuilder.cpp
  )

#
# Adaptors files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Adapters)

APPEND_SET(HEADERS
  Adapters/TSFAnasaziAdapter.hpp
  Adapters/TSF_NVector.hpp
  )

APPEND_SET(SOURCES
  Adapters/TSF_NVector.cpp
  )

#
# Preconditioners files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Preconditioners)

APPEND_SET(HEADERS
  Preconditioners/TSFGenericLeftPreconditioner.hpp
  Preconditioners/TSFGenericRightPreconditioner.hpp
  Preconditioners/TSFILUKPreconditionerFactory.hpp
  Preconditioners/TSFParameterListPreconditionerFactory.hpp
  Preconditioners/TSFPreconditioner.hpp
  Preconditioners/TSFPreconditionerBase.hpp
  Preconditioners/TSFPreconditionerFactory.hpp
  Preconditioners/TSFPreconditionerFactoryBase.hpp 
  )

APPEND_SET(SOURCES
  Preconditioners/TSFParameterListPreconditionerFactory.cpp
  )



IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  TSFExtended
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS sundanceUtils
  )
