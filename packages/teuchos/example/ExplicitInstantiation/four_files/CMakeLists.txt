
INCLUDE(PackageAddTest)

SET(COMMON_SOURCES evilMain.cpp EvilBase.cpp AEvil.cpp BEvil.cpp)

#
# EvilFourFileImplicit
#

PACKAGE_ADD_EXECUTABLE(
  EvilFourFileImplicit
  SOURCES ${COMMON_SOURCES}
  COMM serial mpi
  )

PACKAGE_ADD_TEST(
  EvilFourFileImplicit
  COMM serial mpi
  #PASS_REGULAR_EXPRESSION "blobby"
  )

#
# EvilFourFileExplicit
#

PACKAGE_ADD_EXECUTABLE(
  EvilFourFileExplicit
  SOURCES ${COMMON_SOURCES}
  COMM serial mpi
  )
IF (TARGET Teuchos_EvilFourFileExplicit)
  SET_TARGET_PROPERTIES(
    Teuchos_EvilFourFileExplicit
    PROPERTIES
      COMPILE_DEFINITIONS "DO_EXPLICIT_INSTANTIATION"
    )
ENDIF()

PACKAGE_ADD_TEST(
  EvilFourFileExplicit
  COMM serial mpi
  )

# 2008/08/02: rabartl: TODO: The above calls to PACKAGE_ADD_TEST(...) 
# do not find the test executables and I don't know why.
