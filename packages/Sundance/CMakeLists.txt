
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Sundance  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#



ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Python
  HAVE_SUNDANCE_PYTHON
  "Enable python functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MOOCHO
  HAVE_SUNDANCE_MOOCHO
  "Enable Sundance/Moocho."
  ON )


ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Stokhos
  HAVE_SUNDANCE_STOKHOS
  "Enable Sundance/Stokhos."
  ON )


ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_BROKEN_CODE
  HAVE_SUNDANCE_BROKEN_CODE
  "Enable broken code. For expert developer use only."
  OFF )


ADD_OPTION_AND_DEFINE(TPL_ENABLE_ExodusII
  HAVE_SUNDANCE_EXODUS
  "Enable ExodusII reader in Sundance (requires ExodusII TPL)."
  OFF )

ADD_OPTION_AND_DEFINE(Sundance_ENABLE_EXTERNAL_CHACO
  HAVE_SUNDANCE_CHACO
  "Enable file-based chaco partitioning in Sundance (requires chaco executable in path)."
  OFF )

# ========== Peano configuration ===============
IF(TPL_ENABLE_Peano)
  ADD_OPTION_AND_DEFINE(TPL_ENABLE_Peano
     HAVE_SUNDANCE_PEANO
    "Enable Peano Mesh in Sundance (requires Peano TPL)."
    OFF )  
  INCLUDE_DIRECTORIES(${Peano_INCLUDE_DIRS}) 
  
  IF (TPL_ENABLE_Peano_NO_2D)
    ADD_OPTION_AND_DEFINE(TPL_ENABLE_Peano_NO_2D
      HAVE_SUNDANCE_PEANO_NO_2D
     "Turn OFF Peano 2D"
     ON )
  ELSE()
    # LINK_DIRECTORIES( ${Peano_LIBRARY_DIRS}/libpeano2D.a )
    # PEANO_ADD_LIBRARY(peano2D)
  ENDIF()
  
  IF (TPL_ENABLE_Peano_NO_3D)
    ADD_OPTION_AND_DEFINE(TPL_ENABLE_Peano_NO_3D
      HAVE_SUNDANCE_PEANO_NO_3D
     "Turn OFF Peano 3D"
     ON )
  ELSE()
    # LINK_DIRECTORIES( ${Peano_LIBRARY_DIRS}/libpeano3D.a )
    # PEANO_ADD_LIBRARY(peano3D) 
  ENDIF()
     
ENDIF()
# ============ end of Peano configuration ================

IF(${Proj}_VERBOSE_CONFIGURE)
  PRINT_VAR(HAVE_SUNDANCE_EXODUS)
  PRINT_VAR(HAVE_SUNDANCE_PYTHON)
  PRINT_VAR(HAVE_SUNDANCE_CHACO)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src-utils)
PACKAGE_ADD_TEST_DIRECTORIES(tests-utils)

ADD_SUBDIRECTORY(src-solvers)
PACKAGE_ADD_TEST_DIRECTORIES(tests-solvers)

ADD_SUBDIRECTORY(src-interop)

ADD_SUBDIRECTORY(src-core)
PACKAGE_ADD_TEST_DIRECTORIES(tests-core)

ADD_SUBDIRECTORY(src-std-mesh)
PACKAGE_ADD_TEST_DIRECTORIES(tests-std-mesh)

ADD_SUBDIRECTORY(src-std-framework)
PACKAGE_ADD_TEST_DIRECTORIES(tests-std-framework)

IF(Sundance_ENABLE_Python)
  ADD_SUBDIRECTORY(python)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(examples-tutorial)


#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
