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
#include "itkActiveShapeModel.h"
#include "itkASMPointSampler.h"
#include "itkPointsLocator.h"
#include <vector>
#include "itkStandardMeshRepresenter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTriangleMeshAdapter.h"
#include "itkTransformMeshFilter.h"
#include "Sampling/Logger/BestMatchLogger.h"
#include "MeshOperations.h"
#include "MHFitting.h"
#include "Types.h"

namespace itk {

    typedef itk::Image<float, 3> ImageType;

    typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
    typedef itk::VersorRigid3DTransform<float> RigidTransformType;

    class itkMeshOperations : public mhfitting::MeshOperations<RepresenterType::DatasetPointerType, RepresenterType::PointType> {
        typedef itk::PointsLocator< typename RepresenterType::MeshType::PointsContainer > PointsLocatorType;
        typedef TriangleMeshAdapter<typename RepresenterType::MeshType::PixelType> MeshAdapterType;
        typedef MeshAdapterType::PointNormalType PointNormalType;
        typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
        typedef RepresenterType::PointType PointType;
        typedef statismo::StatisticalModel<RepresenterType::DatasetType>  StatisticalModelType;


    public:
        itkMeshOperations(const StatisticalModelType* statisticalModel, ImageType::Pointer image, PointType rotationCenter) : m_image(image), m_interpolator(InterpolatorType::New()) {
            m_rotationCenter = rotationCenter;
            m_interpolator->SetInputImage(image);
            m_statisticalModel = statisticalModel;
        }


        virtual mhfitting::MHFittingParameters rigidICP(const mhfitting::MHFittingParameters& fittingParameters, std::vector<PointType>) const {
                return fittingParameters;
        }

        void InitializePointLocator(RepresenterType::MeshPointerType mesh) const {
            m_ptLocator = PointsLocatorType::New();
            m_ptLocator->SetPoints(mesh->GetPoints());
            m_ptLocator->Initialize();
        };

        virtual std::pair<RepresenterType::PointType, long> findClosestPoint(RepresenterType::MeshPointerType mesh, RepresenterType::PointType pt) const {
            if (mesh.GetPointer() != m_lastUsedMeshAddress) {
                m_lastUsedMeshAddress = mesh.GetPointer();
                InitializePointLocator(mesh);
            }

            long ptId = m_ptLocator->FindClosestPoint(pt);
            return std::make_pair(mesh->GetPoints()->GetElement(ptId), ptId);
        }

        virtual statismo::VectorType normalAtPoint(RepresenterType::MeshPointerType mesh, RepresenterType::PointType pt) const {
            if (m_normals.GetPointer() == 0) {
                MeshAdapterType::Pointer meshAdapter = MeshAdapterType::New();
                meshAdapter->SetMesh(mesh);
                m_normals = meshAdapter->GetPointNormals();
            }

            long id = this->findClosestPoint(mesh, pt).second;
            PointNormalType normal = m_normals->GetElement(id);
            statismo::VectorType normalAsStatismoVec(3) ;
            normalAsStatismoVec(0) =  normal.GetElement(0);//
            normalAsStatismoVec(1) =  normal.GetElement(1);//
            normalAsStatismoVec(2) =  normal.GetElement(2);//
            return normalAsStatismoVec;
        }


        virtual short huAtPoint(RepresenterType::MeshPointerType mesh, long id) const {
            RepresenterType::PointType pt = mesh->GetPoint(id);
            if (m_interpolator->IsInsideBuffer(pt))
                return m_interpolator->Evaluate(pt);
            else
                return 0;

        }

        RigidTransformType::Pointer getTransform(const mhfitting::MHFittingParameters& fittingParameters) const {

            RigidTransformType::Pointer newRigidTransform = RigidTransformType::New();
            newRigidTransform->SetCenter(m_rotationCenter);
            statismo::VectorType rigidParams = fittingParameters.GetRigidTransformParameters();
            RigidTransformType::ParametersType p(rigidParams.size());
            for (unsigned i = 0 ; i < rigidParams.size(); ++i) {
                p.SetElement(i, rigidParams[i]);
            }
            newRigidTransform->SetParameters(p);
            return newRigidTransform;
        }

        RepresenterType::DatasetPointerType transformMesh(const mhfitting::MHFittingParameters& fittingParameters) const {

            typedef  RepresenterType::DatasetType MeshType ;
            MeshType::Pointer   sampleShape = m_statisticalModel->DrawSample(fittingParameters.GetCoefficients());


            RigidTransformType::Pointer newRigidTransform = getTransform(fittingParameters);


            typedef itk::TransformMeshFilter<MeshType, MeshType, RigidTransformType> TransformMeshFilterType;

            TransformMeshFilterType::Pointer tf = TransformMeshFilterType::New();
            tf->SetInput(sampleShape);
            tf->SetTransform(newRigidTransform);
            tf->Update();

            MeshType::Pointer output = tf->GetOutput();
            return output;

        }




        virtual PointType transformToModelSpace(const statismo::VectorType& rigidTransformParameters, PointType pt) const {
            typedef itk::VersorRigid3DTransform<float> RigidTransformType;


            RigidTransformType::Pointer newRigidTransform = RigidTransformType::New();
            newRigidTransform->SetCenter(m_rotationCenter);
            RigidTransformType::ParametersType p(rigidTransformParameters.size());
            for (unsigned i = 0 ; i < rigidTransformParameters.size(); ++i) {
                p.SetElement(i, rigidTransformParameters[i]);
            }
            newRigidTransform->SetParameters(p);

            return newRigidTransform->GetInverseTransform()->TransformPoint(pt);
        }


        virtual PointType getPointWithId(RepresenterType::MeshPointerType mesh, unsigned id) const {
            return mesh->GetPoint(id);
        }


            private:


        typedef vnl_vector<statismo::ScalarType> VectorType;
        VectorType toVnlVector(const statismo::VectorType& v) {
            return VectorType(v.data(), v.rows());
        }


        ImageType::Pointer m_image;
        InterpolatorType::Pointer m_interpolator;
        PointType m_rotationCenter;
        const StatisticalModelType* m_statisticalModel;
        mutable PointsLocatorType::Pointer m_ptLocator;
        mutable RepresenterType::MeshType* m_lastUsedMeshAddress;
        mutable MeshAdapterType::PointNormalsContainerPointer m_normals;

    };


    template<typename TPointSet, typename TImage>
    class MHFittingResult : public Object {
    public:
        typedef MHFittingResult Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        typedef typename TPointSet::PointType PointType;
        typedef itk::VersorRigid3DTransform<float> RigidTransformType;

        itkNewMacro(Self);
        itkTypeMacro(Self, Object);

        typedef mhfitting::MHFittingParameters ImplType;
        typedef vnl_vector<statismo::ScalarType> VectorType;

        MHFittingResult() : m_meshOperations(0) {}


        void SetInternalData(const itkMeshOperations* meshOperations, ImplType statismoResult) {
            m_samplingParameters = statismoResult;
            m_meshOperations = meshOperations;
        }

        bool IsValid() {
            return true;
        }

//        sampling::MarkovChain< statismo::MHFittingParameters<RigidTransformPointerType> >* GetChain() {
//            return m_statismoResult.GetChain();
//        }


        VectorType GetCoefficients() {
            return toVnlVector(m_samplingParameters.GetCoefficients());
        }

        VersorRigid3DTransform<float>::Pointer GetRigidTransformation() {
            typename itk::VersorRigid3DTransform<float>::Pointer newRigidTransform = itk::VersorRigid3DTransform<float>::New();
            newRigidTransform->SetCenter(GetRotationCenter());
            newRigidTransform->SetParameters(GetRigidTransformParameters());
            return newRigidTransform;
        }



        PointType GetRotationCenter() {
            statismo::VectorType rotationCenter =  m_samplingParameters.GetRotationCenter();
            PointType p;
            assert(rotationCenter.size() == 3);
            for (unsigned i = 0 ; i < 3; ++i) {
                p.SetElement(i, rotationCenter[i]);
            }
            return p;
        }

         RigidTransformType::ParametersType GetRigidTransformParameters() {
             statismo::VectorType rigidParameters =  m_samplingParameters.GetRigidTransformParameters();
             RigidTransformType::ParametersType p(rigidParameters.size());
             for (unsigned i = 0 ; i < rigidParameters.size(); ++i) {
                 p.SetElement(i, rigidParameters[i]);
             }
            return p;
        }





        typename TPointSet::Pointer GetMesh() {
            return m_meshOperations->transformMesh(m_samplingParameters);
        }


    private:
        ImplType m_samplingParameters;
        const itkMeshOperations* m_meshOperations;

        VectorType toVnlVector(const statismo::VectorType& v) {
            return VectorType(v.data(), v.rows());

        }
    };




    template<typename TPointSet, typename TImage>
    class MHFittingStepper : public Object {
    public:
        typedef MHFittingStepper Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(Self, Object);

        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelPointerType;
        typedef typename ActiveShapeModel<TPointSet, TImage>::ImplType::RepresenterType::PointType PointType;
        itk::Rigid3DTransform<float>::Pointer RigidTransformPointerType;
        typedef typename TPointSet::Pointer PointSetPointerType;
        typedef itk::Image<float, 3> ImageType;
        typedef typename itk::Image<float, 3>::Pointer ImagePointerType;
        typedef typename statismo::ASMPreprocessedImage<TPointSet> *PreprocessedImagePointerType;
        typedef mhfitting::MHFittingConfiguration ConfigurationType;
        typedef typename ASMPointSampler<TPointSet, TImage>::Pointer SamplerPointerType;
        typedef mhfitting::MHFittingStep<TPointSet, TImage> ImplType;
        typedef typename mhfitting::BasicSampling<RepresenterType::DatasetType> BasicSamplingType;
        typedef MHFittingResult<TPointSet, TImage> ResultType;

       typedef typename mhfitting::CorrespondencePoints CorrespondencePoints;

        MHFittingStepper() :  m_model(0), m_configuration(statismo::ASMFittingConfiguration(0,0,0)), m_meshOperations(0), m_chain(0) { }

        void init(ImagePointerType targetImage,
                  PreprocessedImagePointerType preprocessedTargetImage,
                  const CorrespondencePoints correspondencePoints,
                  const std::vector<PointType>& targetPoints,
                  ModelPointerType model,
                  SamplerPointerType sampler,
                  ConfigurationType configuration,
                  itk::Rigid3DTransform<float>* transform,
                  statismo::VectorType coeffs)
        {
            m_model = model;
            m_sampler = sampler; // need to hold it here, as otherwise it crashes.
            m_configuration = configuration;
            m_meshOperations = new itkMeshOperations(model->GetStatisticalModel()->GetstatismoImplObj(), targetImage, transform->GetCenter()) ;// TODO this is a memory leak - there must be a better way
            m_preprocessedTargetImage = preprocessedTargetImage;
            m_bestMatchLogger = new sampling::BestMatchLogger<mhfitting::MHFittingParameters>();


            // TODO need pass trough all parameters
//            vector<PointType,int> targetPointsWithIndex,
//            ActiveShapeModelType* asmodel,
//            RigidTransformPointerType transform,
//            statismo::VectorType coeffs

            mhfitting::MHFittingParameters initialParameters(coeffs, fromVnlVector(transform->GetParameters()), fromITKPoint(transform->GetCenter()));

          m_chain = BasicSamplingType::buildInitialPoseChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations, correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_bestMatchLogger);

        }

        std::map<unsigned, statismo::MultiVariateNormalDistribution> computePointUncertainty(const CorrespondencePoints correspondencePoints,
                                     const std::vector<PointType> &targetPoints) {
            mhfitting::MHFittingParameters lastParameters;
            m_chain->current(lastParameters);
            std::cout << "last rigid parameters " << lastParameters.GetRigidTransformParameters() << std::endl;
            std::cout << "last sahpe parameters " << lastParameters.GetCoefficients() << std::endl;
            return BasicSamplingType::estimatePointUncertaintyForInitialPoseChain(
                    m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations, correspondencePoints,
                    targetPoints, m_model->GetstatismoImplObj(), lastParameters);

        }

        void SetChainToLmAndHU(const CorrespondencePoints correspondencePoints, const std::vector<PointType>& targetPoints, itk::Rigid3DTransform<float>* transform, statismo::VectorType coeffs) {
            mhfitting::MHFittingParameters initialParameters(coeffs, fromVnlVector(transform->GetParameters()), fromITKPoint(transform->GetCenter()));
            m_chain = BasicSamplingType::buildLmAndHuChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations, correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_bestMatchLogger);
        }


        void SetChainToPoseOnly(const CorrespondencePoints correspondencePoints, const std::vector<PointType>& targetPoints, itk::Rigid3DTransform<float>* transform, statismo::VectorType coeffs) {
            mhfitting::MHFittingParameters initialParameters(coeffs, fromVnlVector(transform->GetParameters()), fromITKPoint(transform->GetCenter()));
            m_chain = BasicSamplingType::buildInitialPoseChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations, correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_bestMatchLogger);
        }

        void SetChainToPoseAndShape(const CorrespondencePoints correspondencePoints, const std::vector<PointType>& targetPoints, itk::Rigid3DTransform<float>* transform, statismo::VectorType coeffs) {
            mhfitting::MHFittingParameters initialParameters(coeffs, fromVnlVector(transform->GetParameters()), fromITKPoint(transform->GetCenter()));
            m_chain = BasicSamplingType::buildPoseAndShapeChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations, correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_bestMatchLogger);
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

            mhfitting::MHFittingParameters currentParameters = impl->Perform();
            mhfitting::MHFittingParameters bestParameters = m_bestMatchLogger->getBest();
            m_result = ResultType::New();
            m_result->SetInternalData(m_meshOperations, bestParameters);

            delete impl;
        }

        typename ResultType::Pointer GetOutput() {
            return m_result;
        }

    private:
        statismo::VectorType fromVnlVector(const vnl_vector<float>& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

        }

        statismo::VectorType fromITKPoint(const PointType& p) {
            statismo::VectorType v(3);
            v(0) = p.GetElement(0); v(1) = p.GetElement(1); v(2) = p.GetElement(2);
            return v;
        }


    private:
        ModelPointerType m_model;
        ConfigurationType m_configuration;
        itkMeshOperations* m_meshOperations;
        sampling::MarkovChain<mhfitting::MHFittingParameters>* m_chain; // FIXME change type
        PreprocessedImagePointerType m_preprocessedTargetImage;
        SamplerPointerType m_sampler;

        sampling::BestMatchLogger<mhfitting::MHFittingParameters>* m_bestMatchLogger;
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
