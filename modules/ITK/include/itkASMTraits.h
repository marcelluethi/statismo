/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Langguth (christoph.langguth@unibas.ch)
 *
 * Copyright (c) 2011-2015 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef STATISMO_ITKASMTRAITS_H
#define STATISMO_ITKASMTRAITS_H

#include "ActiveShapeModel.h"
#include "itkStandardMeshRepresenter.h"
#include "itkASMFeatureExtractor.h"
#include "itkASMPointSampler.h"
#include "itkActiveShapeModel.h"
#include "itkTriangleMeshAdapter.h"
#include "itkASMFitter.h"

namespace itk {

    template<typename MPixelType, typename IPixelType, unsigned int Dimensions>
    class ASMTraits {
    public:
        typedef MPixelType MeshPixelType;
        typedef IPixelType ImagePixelType;

        static unsigned int GetDimensions() {
            return Dimensions;
        }

        typedef itk::ASMTraits<MeshPixelType, ImagePixelType, Dimensions> ASM;

        typedef itk::Point <MeshPixelType, Dimensions> PointType;
        typedef itk::Vector <MeshPixelType, Dimensions> VectorType;


        typedef itk::Mesh <MeshPixelType, Dimensions> MeshType;
        typedef itk::SmartPointer<MeshType> MeshPointerType;

        typedef itk::Image <ImagePixelType, Dimensions> ImageType;
        typedef itk::SmartPointer<ImageType> ImagePointerType;

        typedef itk::ASMPointSampler<ASM> PointSamplerType;
        typedef itk::SmartPointer<PointSamplerType> PointSamplerPointerType;

        typedef itk::StandardMeshRepresenter<MeshPixelType , Dimensions> MeshRepresenterType;
        typedef itk::SmartPointer<MeshRepresenterType> MeshRepresenterPointerType;

        typedef itk::StatisticalModel<MeshType> StatisticalModelType;
        typedef itk::SmartPointer<StatisticalModelType> StatisticalModelPointerType;

        typedef itk::ActiveShapeModel<ASM> ActiveShapeModelType;
        typedef itk::SmartPointer<ActiveShapeModelType> ActiveShapeModelPointerType;

        typedef itk::ASMFeatureExtractor<ASM> FeatureExtractorType;
        typedef itk::SmartPointer<FeatureExtractorType> FeatureExtractorPointerType;

        typedef itk::ASMFeatureExtractorFactory<ASM> FeatureExtractorFactoryType;
        typedef itk::SmartPointer<FeatureExtractorFactoryType> FeatureExtractorFactoryPointerType;

        typedef itk::TriangleMeshAdapter<MeshPixelType, Dimensions> MeshAdapterType;
        typedef itk::SmartPointer<MeshAdapterType> MeshAdapterPointerType;

        typedef itk::ASMFitterConfiguration FitterConfigurationType;
        typedef itk::SmartPointer<FitterConfigurationType> FitterConfigurationPointerType;

        typedef itk::ASMFitterResult<ASM> FitterResultType;
        typedef itk::SmartPointer<FitterResultType> FitterResultPointerType;

        typedef itk::ASMFitter<ASM> FitterType;
        typedef itk::SmartPointer<FitterType> FitterPointerType;

    };
}

#endif //STATISMO_ITKASMTRAITS_H
