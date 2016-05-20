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

#ifndef STATISMO_ITKMHFITTING_H
#define STATISMO_ITKMHFITTING_H

#include <itkObject.h>
#include <itkMacro.h>
#include "itkASMFitting.h"
#include "ASMFitting.h"
#include "MHFitting.h"
#include "itkActiveShapeModel.h"
#include "itkASMPointSampler.h"
#include "itkPointsLocator.h"
#include <vector>
#include "itkStandardMeshRepresenter.h"
#include "itkASMNormalDirectionImageValueExtractor.h"
#include "itkASMLinearInterpolatedPreprocessedImage.h"

namespace itk {

    typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;

    class itkMeshClosestPoint : public statismo::ClosestPoint<RepresenterType::DatasetPointerType, RepresenterType::PointType> {
        typedef itk::PointsLocator< typename RepresenterType::MeshType::PointsContainer > PointsLocatorType;

    public:
        itkMeshClosestPoint() {}


        void InitializePointLocator(RepresenterType::MeshPointerType mesh) const {
            m_ptLocator = PointsLocatorType::New();
            m_ptLocator->SetPoints(mesh->GetPoints());
            m_ptLocator->Initialize();
        };

        virtual RepresenterType::PointType findClosestPoint(RepresenterType::MeshPointerType mesh, RepresenterType::PointType pt) const {
            if (mesh.GetPointer() != m_cachedMesh.GetPointer()) {
                m_cachedMesh = mesh;
                InitializePointLocator(mesh);
            }

            long ptId = m_ptLocator->FindClosestPoint(pt);
            return mesh->GetPoints()->GetElement(ptId);
        }

    private:
        mutable PointsLocatorType::Pointer m_ptLocator;
        mutable RepresenterType::MeshPointerType m_cachedMesh;
    };


    template<typename TPointSet, typename TImage>
    class MHFittingResult : public Object {
    public:
        typedef MHFittingResult Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelType;
        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelPointerType;
        typedef typename ActiveShapeModel<TPointSet, TImage>::ImplType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;

        itkNewMacro(Self);
        itkTypeMacro(Self, Object);

        typedef statismo::MHFittingResult<RigidTransformPointerType> ImplType;
        typedef vnl_vector<statismo::ScalarType> VectorType;

        void SetInternalData(ImplType statismoResult, ModelPointerType model) {
            m_model = model;
            m_statismoResult = statismoResult;
        }

        bool IsValid() {
            return true;
        }

        sampling::MarkovChain< statismo::MHFittingResult<RigidTransformPointerType> >* GetChain() {
            return m_statismoResult.GetChain();
        }

        VectorType GetCoefficients() {
            return toVnlVector(m_statismoResult.GetCoefficients());
        }

        RigidTransformPointerType GetRigidTransformation() {
            return m_statismoResult.GetRigidTransform();
        }

        typename TPointSet::Pointer GetMesh() {
            typename TPointSet::Pointer instance = m_model->GetStatisticalModel()->DrawSample(GetCoefficients());
            return m_model->GetstatismoImplObj()->GetRepresenter()->TransformMesh(instance, GetRigidTransformation());
        }


    private:
        ImplType m_statismoResult;
        ModelPointerType m_model;

        VectorType toVnlVector(const statismo::VectorType& v) {
            return VectorType(v.data(), v.rows());

        }
    };




    template<typename TPointSet, typename TImage>
    class MHFittingStep : public Object {
    public:
        typedef MHFittingStep Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(Self, Object);

        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelPointerType;
        typedef typename ActiveShapeModel<TPointSet, TImage>::ImplType::RepresenterType::PointType PointType;
        typedef typename ActiveShapeModel<TPointSet, TImage>::ImplType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
        typedef ASMNormalDirectionImageValueExtractor<TPointSet, TImage> ImageValueExtractorType;
        typedef ASMNormalDirectionImageValueExtractorFactory<TPointSet, TImage> ImageValueExtractorFactoryType;
        typedef ASMLinearInterpolatedPreprocessedImage<TPointSet, TImage> LinearInterpolatedPreprocessedImageType;
        typedef typename TPointSet::Pointer PointSetPointerType;
        typedef TImage ImageType;
        typedef typename ImageType::Pointer ImagePointer;
        typedef typename statismo::ASMPreprocessedImage<TPointSet> *PreprocessedImagePointerType;
        typedef statismo::MHFittingConfiguration ConfigurationType;
        typedef typename ASMPointSampler<TPointSet, TImage>::Pointer SamplerPointerType;
        typedef statismo::MHFittingStep<TPointSet, TImage> ImplType;
        typedef MHFittingResult<TPointSet, TImage> ResultType;

        MHFittingStep() : /* m_model(0), m_target(0), m_configuration(statismo::ASMFittingConfiguration(0,0,0)), m_transform(0) */ m_chain(0) { }

        void init(PreprocessedImagePointerType preprocessedFeatureImage,
                  ImagePointer originalImage,
                  const std::vector<PointType>& targetPoints,
                  ModelPointerType model,
                  SamplerPointerType sampler,
                  ConfigurationType configuration,
                  RigidTransformPointerType transform,
                  statismo::VectorType coeffs)
        {
            m_model = model;
//            m_target = target;
            m_sampler = sampler; // need to hold it here, as otherwise it crashes.
//            m_transform = transform;
//            m_coeffs = coeffs;
//            m_configuration = configuration;


            PreprocessedImagePointerType interpolatedImage = LinearInterpolatedPreprocessedImageType::Create(originalImage);

            auto feFactory = ImageValueExtractorFactoryType::GetInstance();
            auto imageValueExtractor = feFactory->Instantiate(19, 0.5);
           
            // TODO need pass trough all parameters
//            vector<PointType,int> targetPointsWithIndex,
//            ActiveShapeModelType* asmodel,
//            RigidTransformPointerType transform,
//            statismo::VectorType coeffs

            typedef typename statismo::mcmc<TPointSet,TImage>::template BasicSampling<TPointSet> BasicSampingType;
            itkMeshClosestPoint* closestPointEv = new itkMeshClosestPoint() ;// TODO this is a memory leak - there must be a better way
            m_chain = BasicSampingType::buildChain(m_model->GetStatisticalModel()->GetRepresenter(), closestPointEv, configuration, preprocessedFeatureImage, interpolatedImage, imageValueExtractor, targetPoints, model->GetstatismoImplObj(), m_sampler, transform, coeffs);
        }


//        void SetSampler(SamplerPointerType sampler) {
//            m_sampler = sampler;
//        }


//        void SetModel(ModelPointerType model) {
//            m_model = model;
//        }

//        void SetCoefficients(statismo::VectorType coeffs) {
//            m_coeffs = coeffs;
//        }

//        void SetLineConstraints(const std::vector<PointType>& linePoints) {
//            m_linePoints = linePoints;
//        }

//        void SetRigidTransformation(RigidTransformPointerType transform) {
//            m_transform = transform;
//        }
//
//        void SetTarget(ImagePointerType target) {
//            m_target = target;
//        }
//
//        void SetSampler(SamplerPointerType sampler) {
//            m_sampler = sampler;
//        }
//
//        void SetConfiguration(const ConfigurationType &configuration) {
//            m_configuration = configuration;
//        }

        void NextSample() {
            ImplType *impl = ImplType::Create(m_chain);

            statismo::MHFittingResult<RigidTransformPointerType> result = impl->Perform();
            m_result = ResultType::New();
            m_result->SetInternalData(result, m_model);

            delete impl;
        }

        typename ResultType::Pointer GetOutput() {
            return m_result;
        }

    private:
        ModelPointerType m_model;
        sampling::MarkovChain<statismo::MHFittingResult<RigidTransformPointerType> >* m_chain; // FIXME change type
 //       statismo::VectorType m_coeffs;
 //       RigidTransformPointerType m_transform;
  //      ImagePointerType m_target;
        SamplerPointerType m_sampler;
    //    ConfigurationType m_configuration;
    //    std::vector<PointType> m_linePoints;
        typename ResultType::Pointer m_result;
    };

    template<typename TPointSet, typename TImage>
    class MHFitting : public Object {
    public:
        typedef MHFitting Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(Self, Object);
    };
}
#endif //STATISMO_ITKASMFITTING_H

