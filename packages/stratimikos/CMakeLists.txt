
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

PACKAGE(Stratimikos)

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

ADD_SUBDIRECTORY(adapters)

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_POSTPROCESS()

#
# ToDos:
#
# (*) 2008/12/10: rabartl: Add support for Belos once Belos has been
# CMakeifyied
#
# (*) 2008/12/10: rabartl: Add general support for generating sample
# output used by doxygen documentation.
#
# (*) 2008/12/10: rabartl: Change PACKAGE_ADD_LIBRARY(...) to
# accept DEPPACKAGES.  Then, use this with the adatper and main
# libraries in Stratimikos to better define the library dependency
# structure.  Right now, for example, the thyra/amesos adatpers depend
# on ML, AztecOO, etc when they should not.
#
