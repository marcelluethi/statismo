project(RepresenterTests)
cmake_minimum_required(VERSION 2.6)


enable_testing()
#INCLUDE(Dart)

find_package(ITK 4.1 REQUIRED)
include (${ITK_USE_FILE})

find_package(VTK REQUIRED)
include (${VTK_USE_FILE})

string(REGEX REPLACE "itkhdf5_cpp" "" ITK_LIBRARIES "${ITK_LIBRARIES}")
string(REGEX REPLACE "itkhdf5" "" ITK_LIBRARIES "${ITK_LIBRARIES}")

include_directories(
  ${STATISMO_ROOT_DIR}
  ${STATISMO_ROOT_DIR}/statismo_ITK
  ${STATISMO_INCLUDE_DIR}
  ${ITK_INCLUDE_DIRS}
  ${STATISMO_ROOT_DIR}/Representers/ITK
  ${STATISMO_ROOT_DIR}/Representers/VTK
  ${VTK_INCLUDE_DIR}
)

add_executable(itkVectorImageRepresenterTest itkVectorImageRepresenterTest.cpp)
target_link_libraries(itkVectorImageRepresenterTest ${ITK_LIBRARIES} ${STATISMO_LIBRARIES})

add_executable(itkImageRepresenterTest itkImageRepresenterTest.cpp)
target_link_libraries(itkImageRepresenterTest   ${ITK_LIBRARIES} ${STATISMO_LIBRARIES})

add_executable(itkMeshRepresenterTest itkMeshRepresenterTest.cpp)
target_link_libraries(itkMeshRepresenterTest   ${ITK_LIBRARIES} ${STATISMO_LIBRARIES})

add_executable(vtkPolyDataRepresenterTest vtkPolyDataRepresenterTest.cpp)
target_link_libraries(vtkPolyDataRepresenterTest ${VTK_LIBRARIES} ${STATISMO_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(vtkStructuredPointsRepresenterTest vtkStructuredPointsRepresenterTest.cpp)
target_link_libraries(vtkStructuredPointsRepresenterTest ${VTK_LIBRARIES} ${STATISMO_LIBRARIES}  ${HDF5_LIBRARIES})

add_test(itkVectorImageRepresenterTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkVectorImageRepresenterTest ${STATISMO_ROOT_DIR}/data)
add_test(itkImageRepresenterTest  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkImageRepresenterTest ${STATISMO_ROOT_DIR}/data)
add_test(itkMeshRepresenterTest  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkMeshRepresenterTest ${STATISMO_ROOT_DIR}/data)
add_test(vtkPolyDataRepresenterTest  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/vtkPolyDataRepresenterTest ${STATISMO_ROOT_DIR}/data)
add_test(vtkStructuredPointsRepresenterTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/vtkStructuredPointsRepresenterTest ${STATISMO_ROOT_DIR}/data)
