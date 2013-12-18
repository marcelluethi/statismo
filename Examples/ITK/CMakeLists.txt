cmake_minimum_required(VERSION 2.6)
project(ssm_itk)

find_package(ITK 4.1 REQUIRED)
include (${ITK_USE_FILE})

set (ITK_STATISTICAL_MODEL_INCLUDE_DIR ${STATISMO_ROOT_DIR}/ITKstatismo)

file(GLOB INCLUDE_FILES *.h)
file(GLOB SOURCE_FILES  *.txx )

include_directories(
  ${STATISMO_ROOT_DIR}
  ${STATISMO_ROOT_DIR}/statismo_ITK
  ${STATISMO_INCLUDE_DIR}
	${ITK_INCLUDE_DIRS}
  ${STATISMO_ROOT_DIR}/Representers/ITK
)

# link_directories( ${HDF5_LIBRARY_DIR} ${ITK_LIBRARY_DIRS} ${STATISMO_LIBRARY_DIR} )

add_executable( itkBuildDeformationModel
  buildDeformationModel.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkBuildDeformationModel ${ITK_LIBRARIES})

add_executable (itkDeformationModelFitting
  deformationModelFitting.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkDeformationModelFitting   ${ITK_LIBRARIES})

add_executable(itkBuildShapeModel
  buildShapeModel.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkBuildShapeModel ${ITK_LIBRARIES} )

add_executable (itkBuildShapeModel_75pcvar buildShapeModel_75pcvar.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkBuildShapeModel_75pcvar ${ITK_LIBRARIES} )

add_executable (itkBuildImageIntensityModelOnROI buildImageIntensityModelOnROI.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkBuildImageIntensityModelOnROI ${ITK_LIBRARIES} )

add_executable (itkShapeModelFitting shapeModelFitting.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkShapeModelFitting  ${ITK_LIBRARIES})

add_executable (itkLandmarkConstrainedShapeModelFitting landmarkConstrainedShapeModelFitting.cpp
	${SOURCE_FILES}
	${INCLUDE_FILES}
)
target_link_libraries (itkLandmarkConstrainedShapeModelFitting  ${ITK_LIBRARIES})

add_executable (itkBuildGaussianProcessDeformationModel buildGaussianProcessDeformationModel.cpp
        ${SOURCE_FILES}
        ${INCLUDE_FILES}
)
target_link_libraries (itkBuildGaussianProcessDeformationModel ${ITK_LIBRARIES} )

add_executable(itkLowRankGaussianProcessImageToImageRegistration lowRankGaussianProcessImageToImageRegistration.cpp
        ${SOURCE_FILES}
        ${INCLUDE_FILES}
)
target_link_libraries (itkLowRankGaussianProcessImageToImageRegistration  ${ITK_LIBRARIES})

add_executable(itkSimpleGaussianProcessImageToImageRegistration simpleGaussianProcessImageToImageRegistration.cpp
        ${SOURCE_FILES}
        ${INCLUDE_FILES}
)
target_link_libraries(itkSimpleGaussianProcessImageToImageRegistration  ${ITK_LIBRARIES})

INSTALL(TARGETS
itkBuildDeformationModel
itkDeformationModelFitting
itkBuildShapeModel
itkBuildShapeModel_75pcvar
itkShapeModelFitting
itkLandmarkConstrainedShapeModelFitting
itkBuildImageIntensityModelOnROI
itkBuildGaussianProcessDeformationModel
itkLowRankGaussianProcessImageToImageRegistration
itkSimpleGaussianProcessImageToImageRegistration
RUNTIME DESTINATION bin)
INSTALL(PROGRAMS runITKExamples.sh DESTINATION  bin)
