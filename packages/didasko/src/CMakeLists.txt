
#This file is only here to set up the config file in the right spot.
#potentially we could put the latex building logic in here as well at
#some point.

INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

IF(HAVE_DIDASKO_EPETRAEXT AND ${PACKAGE_NAME}_ENABLE_Zoltan)
  SET(HAVE_ZOLTAN ON)
ENDIF()

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h ${PACKAGE_NAME}_ConfigDefs.h
  DESTINATION "include"
)
