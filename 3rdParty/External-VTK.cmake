message( "External project - VTK" )

find_package(Git)
if(NOT GIT_FOUND)
  message(ERROR "Cannot find git. git is required for Superbuild")
endif()

option( USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

set( VTK_DEPENDENCIES HDF5 )

ExternalProject_Add(VTK
  DEPENDS ${VTK_DEPENDENCIES}
  GIT_REPOSITORY ${git_protocol}://vtk.org/VTK.git
  GIT_TAG v6.0.0
  SOURCE_DIR VTK
  BINARY_DIR VTK-build
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DBUILD_TESTING:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DVTK_BUILD_ALL_MODULES:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DITK_USE_SYSTEM_HDF5:BOOL=ON
    -DHDF5_DIR:PATH=${HDF5_DIR}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/3rdParty/VTK
)

set( VTK_DIR ${CMAKE_BINARY_DIR}/3rdParty/VTK/lib/cmake/vtk-6.0/ )
