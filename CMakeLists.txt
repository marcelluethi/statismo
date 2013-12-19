cmake_minimum_required(VERSION 2.8.5)
project(statismo)

enable_testing()

set(STATISMO_VERSION 0.8)
set(STATISMO_MINOR_VERSION 1)


# Setup build locations.
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

set( CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

find_package( Boost REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )


  find_package( HDF5 COMPONENTS C CXX )
  message( "HDF5_LIBRARIES: ${HDF5_LIBRARIES}" )
  message( "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}" )

  include_directories( ${HDF5_INCLUDE_DIRS} )

# rpath settings
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${STATISMO_LIBRARY_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#
#optional examples and wrapping
option(BUILD_ITK_EXAMPLES "Build the itk examples " OFF)
option(BUILD_VTK_EXAMPLES "Build the vtk examples " OFF)
option(BUILD_PYTHON_WRAPPING_ITK "Build python wrapping for the itk representers using wrapitk (experimental)" OFF)
option(BUILD_PYTHON_WRAPPING_VTK "Build python wrapping for the VTK representers" OFF)

mark_as_advanced(BUILD_PYTHON_WRAPPING_ITK)

#
#optional tests
#
option(BUILD_STATISMO_TESTS "Build tests" OFF)
option(BUILD_STATISMO_TESTS_PYTHON "Build Unit tests written in python (requires python wrapping)")
option(BUILD_REPRESENTER_TESTS "Build representer tests (requires ITK and VTK)" OFF)
mark_as_advanced(BUILD_REPRESENTER_TESTS)

option(HAS_CXX11_ASYNC "Statismo can make use of multi-threading in case the compiler supports std::async of the C++11 standard " OFF)
if(HAS_CXX11_ASYNC )
  add_definitions( -DHAS_CXX11_ASYNC )
endif(HAS_CXX11_ASYNC)


#
# Install boost and eigen, by just copying them from the 3rdParty directory
#
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/boost DESTINATION include/)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/Eigen DESTINATION include/)

#
# Install the statismo and statismo itk files (these are header only, so no compilation needed)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/statismo DESTINATION include)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/statismo_ITK DESTINATION include)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Representers DESTINATION include)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/statismo-config.cmake.in  ${CMAKE_CURRENT_BINARY_DIR}/statismo-config.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/statismo-config.cmake DESTINATION lib/cmake/statismo-${STATISMO_VERSION}${STATISMO_MINOR_VERSION})

#
# copy the data
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION share)

set(STATISMO_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_ITK_EXAMPLES)
  add_subdirectory(Examples/ITK)
endif()

if(BUILD_VTK_EXAMPLES)
  add_subdirectory(Examples/VTK)
endif()

if(BUILD_STATISMO_TESTS)
  add_subdirectory(Tests)
endif()

if(BUILD_REPRESENTER_TESTS)
add_subdirectory(Tests/representerTests)
endif()

if(BUILD_STATISMO_TESTS_PYTHON)
  set(BUILD_PYTHON_WRAPPING_VTK ON)
  add_custom_target(unit_test  ${CMAKE_CURRENT_SOURCE_DIR}/Tests/statismoTests/runPythonTests.sh ${CMAKE_CURRENT_BINARY_DIR}/Wrapping/SWIG ${CMAKE_CURRENT_SOURCE_DIR}/Tests/statismoTests/)
endif()


if(BUILD_PYTHON_WRAPPING_VTK)
  add_subdirectory(Wrapping/SWIG)
endif()

if(BUILD_PYTHON_WRAPPING_ITK)
  add_subdirectory(Wrapping/ITK)
endif()
