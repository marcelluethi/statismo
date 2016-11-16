//
// Created by luetma00 on 15.11.16.
//

#ifndef STATISMO_TYPES_H_H
#define STATISMO_TYPES_H_H

#include <itkPoint.h>

namespace mhfitting {

    typedef itk::Image<float, 3> ImageType;
typedef  itk::Point<float, 3> PointType;
    typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
    typedef RepresenterType::DatasetType MeshType;
    typedef RepresenterType::PointType PointType;

    typedef statismo::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
    typedef statismo::ASMPointSampler<MeshType, ImageType> PointSamplerType;
    typedef statismo::StatisticalModel<MeshType> StatisticalModelType;
    typedef statismo::ASMFeatureExtractor<MeshType, ImageType> FeatureExtractorType;
    typedef statismo::ASMPreprocessedImage<MeshType> PreprocessedImageType;
    typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
    typedef  statismo::ASMFittingStep<MeshType, ImageType> ASMFittingStepType;
    typedef std::vector<std::pair<unsigned, PointType> > CorrespondencePoints;

}

#endif //STATISMO_TYPES_H_H
