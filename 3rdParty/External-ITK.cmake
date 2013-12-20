message( "External project - ITK" )

find_package(Git)
if(NOT GIT_FOUND)
  message(ERROR "Cannot find git. git is required for Superbuild")
endif()

option( USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

set( ITK_DEPENDENCIES HDF5 VTK )

ExternalProject_Add(ITK
  DEPENDS ${ITK_DEPENDENCIES}
  GIT_REPOSITORY ${git_protocol}://itk.org/ITK.git
  GIT_TAG v4.5.0
  SOURCE_DIR ITK
  BINARY_DIR ITK-build
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DBUILD_TESTING:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
    -DModule_ITKVtkGlue:BOOL=ON
    -DModule_ITKReview:BOOL=ON
    -DITK_LEGACY_REMOVE:BOOL=ON
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/3rdParty/ITK
    -DITK_USE_SYSTEM_HDF5:BOOL=ON
    -DHDF5_DIR:PATH=${HDF5_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
)

set( ITK_DIR ${CMAKE_BINARY_DIR}/3rdParty/ITK/lib/cmake/ITK-4.5/ )
