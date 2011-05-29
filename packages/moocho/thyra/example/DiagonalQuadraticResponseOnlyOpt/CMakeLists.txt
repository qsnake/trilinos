INCLUDE(PackageAddExecutableAndTest)


#
# EpetraExtDiagonalQuadraticResponseOnlyOpt
#


SET(EPETRAEXT_DQROO_DIR
  ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalQuadraticResponseOnlyOpt)


INCLUDE_DIRECTORIES(${EPETRAEXT_DQROO_DIR})


PACKAGE_ADD_EXECUTABLE_AND_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  SOURCES
    EpetraExtDiagonalQuadraticResponseOnlyOptMain.cpp
    ${EPETRAEXT_DQROO_DIR}/EpetraExt_DiagonalQuadraticResponseOnlyModelEvaluator.cpp
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np1
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np2
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  COMM mpi
  NUM_MPI_PROCS 2
  )


PACKAGE_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np3
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  COMM mpi
  NUM_MPI_PROCS 3
  )


PACKAGE_ADD_TEST( EpetraExtDiagonalQuadraticResponseOnlyOpt
  NAME_POSTFIX np4
  ARGS "--local-dim=100 --pt=0.0 --p0=2.0"
  COMM mpi
  NUM_MPI_PROCS 4
  )


#
# ThyraDiagonalQuadraticResponseOnlyOpt
#


PACKAGE_ADD_EXECUTABLE_AND_TEST( ThyraDiagonalQuadraticResponseOnlyOpt
  SOURCES
    ThyraDiagonalQuadraticResponseOnlyOptMain.cpp
  DEPLIBS thyra_nonlinear_test_models
  NUM_MPI_PROCS 1
  )

# ToDo: Add some more tests ...
