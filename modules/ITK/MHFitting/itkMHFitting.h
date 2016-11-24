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
#include <itkPosteriorModelBuilder.h>
#include "Sampling/Logger/BestMatchLogger.h"
#include "MeshOperations.h"
#include "MHFitting.h"
#include "Types.h"
#include "MHFittingUtils.h"
#include "Sampling/Logger/LoggerDistributor.h"
#include "Sampling/Logger/AcceptRejectLogger.h"




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
        typedef itk::PosteriorModelBuilder<typename RepresenterType::DatasetType> PosteriorModelBuilderType;
        typedef statismo::StatisticalModel<RepresenterType::DatasetType>  StatisticalModelType;


    public:
        itkMeshOperations(const StatisticalModelType* statisticalModel, const mhfitting::CenterLineModel& centerlineModel, ImageType::Pointer image, PointType rotationCenter)
                :   m_image(image), m_interpolator(InterpolatorType::New()), m_centerLineModel(centerlineModel) {
            m_rotationCenter = rotationCenter;
            m_interpolator->SetInputImage(image);
            m_statisticalModel = statisticalModel;
            m_ptLocator = PointsLocatorType::New();



        }


//
//        virtual mhfitting::MHFittingParameters::ModelParameters projectModelInstanceToTargetPoints(const mhfitting::MHFittingParameters &fittingParameters,
//                                                                       const mhfitting::CorrespondencePoints& correspondencePoints,
//                                                                       const std::vector<PointType> &targetPoints) const {
//
//            RepresenterType::MeshType::Pointer centerLine = transformMesh(fittingParameters);
//            RigidTransformType::Pointer rigidTransform = getTransformFromParameters(fittingParameters);
//
//            // transform all the points into model space and add them as constraint
//
//            PosteriorModelBuilderType::PointValueListType constraints;
//
//            for (std::vector<PointType>::const_iterator ptIt = targetPoints.begin(); ptIt != targetPoints.end(); ++ptIt) {
//                long pointId = findClosestPoint(centerLine, *ptIt).second;
//                PointType refPt = m_statisticalModel->GetRepresenter()->GetReference()->GetPoint(pointId);
//                constraints.push_back(std::make_pair(refPt, rigidTransform->GetInverseTransform()->TransformPoint(*ptIt)));
//            }
//
//            for (mhfitting::CorrespondencePoints::const_iterator cpIt = correspondencePoints.begin(); cpIt != correspondencePoints.end(); ++cpIt) {
//                PointType refPt = m_statisticalModel->GetRepresenter()->GetReference()->GetPoint(cpIt->first);
//                constraints.push_back(std::make_pair(refPt, rigidTransform->GetInverseTransform()->TransformPoint(cpIt->second)));
//            }
//
//            statismo::VectorType newCoefficients = m_statisticalModel->ComputeCoefficientsForPointValues(constraints, 2.0);
//
//
//            return mhfitting::MHFittingParameters::ModelParameters(newCoefficients);
//
//        }


        virtual mhfitting::MHFittingParameters::ModelParameters projectModelInstanceToTargetPoints(const mhfitting::MHFittingParameters &fittingParameters,
                                                                                                   const mhfitting::CorrespondencePoints& correspondencePoints,
                                                                                                   const std::vector<PointType> &targetPoints) const {


            RepresenterType::MeshType::Pointer centerLine = transformCenterLine(fittingParameters);
            RigidTransformType::Pointer rigidTransform = getTransformFromParameters(fittingParameters);

            // transform all the points into model space and add them as constraint

            PosteriorModelBuilderType::PointValueListType constraints;
            using mhfitting::MHFittingUtils;

            // compute center
            std::list<PointType> centerPointsTarget;
            std::list<PointType> centerPointsModel;

            MHFittingUtils::LineMapType lineMap = MHFittingUtils::getLineMap(targetPoints);
            for (MHFittingUtils::LineMapType::const_iterator it = lineMap.begin(); it != lineMap.end(); ++it) {
                std::vector<PointType> linePoints = it->second;
                double centerX =0; double centerY = 0; double centerZ = 0;
                for (std::vector<PointType>::const_iterator ptIt = linePoints.begin(); ptIt != linePoints.end(); ++ptIt) {
                    centerX += ptIt->GetElement(0) / linePoints.size();
                    centerY += ptIt->GetElement(1) / linePoints.size();
                    centerZ += ptIt->GetElement(2) / linePoints.size();
                }
                PointType centerPoint;
                centerPoint.SetElement(0, centerX); centerPoint.SetElement(1, centerY); centerPoint.SetElement(2, centerZ);
                centerPointsTarget.push_back(centerPoint);

                long pointId = findClosestPoint(centerLine, centerPoint).second;

                PointType refPt = m_centerLineModel.GetStatisticalModel()->GetRepresenter()->GetReference()->GetPoint(pointId);
                centerPointsModel.push_back(refPt);
                constraints.push_back(std::make_pair(refPt, rigidTransform->GetInverseTransform()->TransformPoint(centerPoint)));
            }



            mhfitting::UIObject& ui = mhfitting::UIObject::getInstance();
      //      ui.ui.showPointCloud(ui.targetGroup, centerPointsTarget, "target center");
    //            ui.ui.showPointCloud(ui.modelGroup, centerPointsModel, "model center");

      //      ui.ui.showTriangleMesh(ui.targetGroup, m_centerLineModel.GetStatisticalModel()->GetRepresenter()->GetReference(), "centerline Ref");
    //            ui.ui.showTriangleMesh(ui.targetGroup, centerLine, "centerline current");
            statismo::VectorType newCoefficients = m_centerLineModel.GetStatisticalModel()->ComputeCoefficientsForPointValues(constraints, 1e-1);

            return mhfitting::MHFittingParameters::ModelParameters(newCoefficients);

        }
        virtual mhfitting::MHFittingParameters::RigidParameters alignModelInstanceToTargetPoints(
                const mhfitting::MHFittingParameters &fittingParameters, const std::vector<PointType> &targetPoints) const {

            PointType center;

            statismo::VectorType rotationCenter = fittingParameters.GetRotationCenter().GetParameters();
            for (unsigned i = 0; i < 3; ++i) {center[i] = rotationCenter[i]; }

            // find the closest points to harget points
            RepresenterType::MeshType::Pointer currentMesh = transformMesh(fittingParameters);
            RepresenterType::MeshType::Pointer nonAlignedMesh = m_statisticalModel->DrawSample(fittingParameters.GetCoefficients().GetParameters());

            std::vector<PointType> closestPointsOnNonAlignedMesh;
            //std::vector<PointType> closestPointOnCurrentMesh;
            for (std::vector<PointType>::const_iterator ptIt = targetPoints.begin(); ptIt != targetPoints.end(); ++ptIt) {
                long ptId = findClosestPoint(currentMesh, *ptIt).second;
                closestPointsOnNonAlignedMesh.push_back(nonAlignedMesh->GetPoint(ptId));
            }


            // create fake rotation
            /*
            statismo::VectorType fakeRotParams = statismo::VectorType::Zero(3); fakeRotParams(0) = 3.14 / 4.0; fakeRotParams(1) = 1.0; fakeRotParams(2) = 3.14 / 8.0;

            mhfitting::MHFittingParameters fakeParams(
                    mhfitting::MHFittingParameters::ModelParameters(statismo::VectorType::Zero(m_statisticalModel->GetNumberOfPrincipalComponents())),
                    mhfitting::MHFittingParameters::RotationParameters(fakeRotParams),
                    mhfitting::MHFittingParameters::TranslationParameters(statismo::VectorType::Zero(3)),
                    mhfitting::MHFittingParameters::RotationCenter(statismo::VectorType::Zero(3))
            );
            VersorRigid3DTransform<float>::Pointer rigidTransformx = getTransformFromParameters(fakeParams);

            std::vector<PointType> fakePoints;
            for (std::vector<PointType>::const_iterator ptIt = targetPoints.begin(); ptIt != targetPoints.end(); ++ptIt) {
                fakePoints.push_back(rigidTransformx->TransformPoint(*ptIt));
            }
            std::cout << " parameters fake transform " << rigidTransformx->GetParameters() << std::endl;
            */
            mhfitting::MHFittingUtils::ProcrustesResult pResult = mhfitting::MHFittingUtils::computeProcrustes(closestPointsOnNonAlignedMesh, targetPoints, center);


            // We need to transform things back into a versor. Therefore we create a versor transform set its matrix
            // and hope for the best that we get a versor back. We need double precision, as otherwise ITK complains.
            VersorRigid3DTransform<double>::Pointer rigidTransform = VersorRigid3DTransform<double>::New();
            VersorRigid3DTransform<double>::MatrixType R;
            VersorRigid3DTransform<double>::VectorType t;

            for (unsigned i = 0; i < 3; ++i) {
                t[i] = pResult.translationVector(i);
                for (unsigned j = 0; j < 3; j++) {
                    R(i,j) = pResult.rotationMatrix(i,j);
                }
            }


            rigidTransform->SetCenter(center);
            rigidTransform->SetMatrix(R);
            rigidTransform->SetTranslation(t);
            rigidTransform->Modified();
            statismo::VectorType newRotationParameters(3);

            for (unsigned i = 0; i < 3; ++i) {
                newRotationParameters(i) = rigidTransform->GetParameters().GetElement(i); // the first three entries are the versor components
            }

            statismo::VectorType newTranslationParameters = pResult.translationVector.cast<float>();
            mhfitting::MHFittingParameters::RigidParameters newParameters(mhfitting::MHFittingParameters::RotationParameters(newRotationParameters), mhfitting::MHFittingParameters::TranslationParameters(newTranslationParameters), fittingParameters.GetRotationCenter());

            return newParameters;
        }

        void InitializePointLocator(RepresenterType::MeshPointerType mesh) const {
            if (m_lastUsedMesh.GetPointer() == 0 || mesh->GetPoints() !=  m_lastUsedMesh->GetPoints()) {
                m_lastUsedMesh = mesh;
                m_ptLocator->SetPoints(mesh->GetPoints());
                m_ptLocator->Initialize();
           }
        };

        virtual std::pair<RepresenterType::PointType, long> findClosestPoint(RepresenterType::MeshPointerType mesh, RepresenterType::PointType pt) const {

            InitializePointLocator(mesh);
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

        RigidTransformType::Pointer getTransformFromParameters(const mhfitting::MHFittingParameters &fittingParameters) const {

            RigidTransformType::Pointer newRigidTransform = RigidTransformType::New();
            newRigidTransform->SetCenter(m_rotationCenter);
            statismo::VectorType rotationParams = fittingParameters.GetRotationParameters().GetParameters();
            statismo::VectorType translationParams = fittingParameters.GetTranslationParameters().GetParameters();
            RigidTransformType::ParametersType p(rotationParams.size() + translationParams.size());
            for (unsigned i = 0 ; i < rotationParams.size(); ++i) {
                p.SetElement(i, rotationParams[i]);
            }
            for (unsigned i = 0 ; i < translationParams.size(); ++i) {
                p.SetElement(i + 3, translationParams[i]);
            }

            newRigidTransform->SetParameters(p);
            return newRigidTransform;
        }


        mhfitting::MHFittingParameters getParametersFromTransform(RigidTransformType *rigidTransform,
                                                                  statismo::VectorType coeffs) const {

            typename RigidTransformType::ParametersType rigidTransformParams = rigidTransform->GetParameters();
            statismo::VectorType translationParams(3);
            statismo::VectorType rotationParams(3);
            statismo::VectorType rotationCenter(3);

            for (unsigned i = 0; i < 3; ++i) {
                rotationParams(i) = rigidTransformParams[i];
                rotationCenter(i) = rigidTransform->GetCenter().GetElement(i);
                translationParams(i) = rigidTransformParams[3 + i];
            }

            return mhfitting::MHFittingParameters(
                    mhfitting::MHFittingParameters::ModelParameters(coeffs),
                    mhfitting::MHFittingParameters::RotationParameters(rotationParams),
                    mhfitting::MHFittingParameters::TranslationParameters(translationParams),
                    mhfitting::MHFittingParameters::RotationCenter(rotationCenter)
            );

            //return params;
        }


        RepresenterType::DatasetPointerType transformMesh(const mhfitting::MHFittingParameters& fittingParameters) const {

            return transformMeshInternal(m_statisticalModel, fittingParameters);
        }


        RepresenterType::DatasetPointerType transformCenterLine(const mhfitting::MHFittingParameters& fittingParameters) const {

          return transformMeshInternal(m_centerLineModel.GetStatisticalModel(), fittingParameters);

        }

        RepresenterType::DatasetPointerType transformMeshInternal(const StatisticalModelType* model, const mhfitting::MHFittingParameters& fittingParameters) const {

            typedef  RepresenterType::DatasetType MeshType ;
            MeshType::Pointer   sampleShape = model->DrawSample(fittingParameters.GetCoefficients().GetParameters());


            RigidTransformType::Pointer newRigidTransform = getTransformFromParameters(fittingParameters);


            typedef itk::TransformMeshFilter<MeshType, MeshType, RigidTransformType> TransformMeshFilterType;

            TransformMeshFilterType::Pointer tf = TransformMeshFilterType::New();
            tf->SetInput(sampleShape);
            tf->SetTransform(newRigidTransform);
            tf->Update();

            MeshType::Pointer output = tf->GetOutput();
            return output;

        }




        virtual PointType transformToModelSpace(const mhfitting::MHFittingParameters& fittingParameters, PointType pt) const {
            typedef itk::VersorRigid3DTransform<float> RigidTransformType;

            RigidTransformType::Pointer rigidTransform = getTransformFromParameters(fittingParameters);

            return rigidTransform->GetInverseTransform()->TransformPoint(pt);
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
        mhfitting::CenterLineModel m_centerLineModel;
        mutable PointsLocatorType::Pointer m_ptLocator;
        mutable RepresenterType::MeshType::Pointer m_lastUsedMesh;
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
            return toVnlVector(m_samplingParameters.GetCoefficients().GetParameters());
        }

        VersorRigid3DTransform<float>::Pointer GetRigidTransformation() {
            typename itk::VersorRigid3DTransform<float>::Pointer newRigidTransform = itk::VersorRigid3DTransform<float>::New();
            newRigidTransform->SetCenter(GetRotationCenter());
            newRigidTransform->SetParameters(GetRigidTransformParameters());
            return newRigidTransform;
        }



        PointType GetRotationCenter() {
            statismo::VectorType rotationCenter =  m_samplingParameters.GetRotationCenter().GetParameters();
            PointType p;
            assert(rotationCenter.size() == 3);
            for (unsigned i = 0 ; i < 3; ++i) {
                p.SetElement(i, rotationCenter[i]);
            }
            return p;
        }

         RigidTransformType::ParametersType GetRigidTransformParameters() {
             return  m_meshOperations->getTransformFromParameters(m_samplingParameters)->GetParameters();
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

        MHFittingStepper() :  m_model(0), m_configuration(statismo::ASMFittingConfiguration(0,0,0)) { }

        void init(ImagePointerType targetImage,
                  PreprocessedImagePointerType preprocessedTargetImage,
                  const CorrespondencePoints correspondencePoints,
                  const std::vector<PointType>& targetPoints,
                  ModelPointerType model,
                  const mhfitting::CenterLineModel& centerLineModel,
                  SamplerPointerType sampler,
                  ConfigurationType configuration,
                  RigidTransformType* transform,
                  statismo::VectorType coeffs,
                  std::ofstream& logStream)
        {
            m_model = model;
            m_sampler = sampler; // need to hold it here, as otherwise it crashes.
            m_configuration = configuration;
            m_meshOperations.reset(new itkMeshOperations(model->GetStatisticalModel()->GetstatismoImplObj(), centerLineModel, targetImage, transform->GetCenter()));
            m_preprocessedTargetImage = preprocessedTargetImage;

            m_bestMatchLogger.reset(new sampling::BestMatchLogger<mhfitting::MHFittingParameters>());
            m_loggerDistributor.reset(new sampling::LoggerDistributor<mhfitting::MHFittingParameters>());
            m_printLogger.reset(new sampling::AcceptRejectLogger<mhfitting::MHFittingParameters>(logStream));
            m_loggerDistributor->add(m_bestMatchLogger.get());
            m_loggerDistributor->add(m_printLogger.get());

            // TODO need pass trough all parameters
//            vector<PointType,int> targetPointsWithIndex,
//            ActiveShapeModelType* asmodel,
//            RigidTransformPointerType transform,
//            statismo::VectorType coeffs
            mhfitting::MHFittingParameters initialParameters = m_meshOperations->getParametersFromTransform(transform, coeffs);

          m_chain.reset(BasicSamplingType::buildInitialPoseChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations.get(), correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_loggerDistributor.get()));

        }

//        std::map<unsigned, statismo::MultiVariateNormalDistribution> computePointUncertainty(const CorrespondencePoints correspondencePoints,
//                                     const std::vector<PointType> &targetPoints) {
//            mhfitting::MHFittingParameters lastParameters;
//            m_chain->current(lastParameters);
//
//            return BasicSamplingType::estimatePointUncertaintyForInitialPoseChain(
//                    m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations, correspondencePoints,
//                    targetPoints, m_model->GetstatismoImplObj(), lastParameters);
//
//        }

        void SetChainToLmAndHU(const CorrespondencePoints correspondencePoints, const std::vector<PointType>& targetPoints, RigidTransformType* transform, statismo::VectorType coeffs) {
            mhfitting::MHFittingParameters initialParameters = m_meshOperations->getParametersFromTransform(transform, coeffs);
            m_chain.reset(BasicSamplingType::buildLmAndHuChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations.get(), correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_loggerDistributor.get()));
        }


        void SetChainToPoseOnly(const CorrespondencePoints correspondencePoints, const std::vector<PointType>& targetPoints, RigidTransformType* transform, statismo::VectorType coeffs) {
            mhfitting::MHFittingParameters initialParameters = m_meshOperations->getParametersFromTransform(transform, coeffs);
            m_chain.reset(BasicSamplingType::buildInitialPoseChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations.get(), correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_loggerDistributor.get()));
        }

        void SetChainToPoseAndShape(const CorrespondencePoints correspondencePoints, const std::vector<PointType>& targetPoints, RigidTransformType* transform, statismo::VectorType coeffs) {
            mhfitting::MHFittingParameters initialParameters = m_meshOperations->getParametersFromTransform(transform, coeffs);
            m_chain.reset(BasicSamplingType::buildPoseAndShapeChain(m_model->GetStatisticalModel()->GetRepresenter(), m_meshOperations.get(), correspondencePoints, targetPoints, m_model->GetstatismoImplObj(), initialParameters, m_loggerDistributor.get()));
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
            ImplType *impl = ImplType::Create(m_chain.get());

            mhfitting::MHFittingParameters currentParameters = impl->Perform();
            mhfitting::MHFittingParameters bestParameters = m_bestMatchLogger->getBest();
            m_result = ResultType::New();
            m_result->SetInternalData(m_meshOperations.get(), bestParameters);

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
        std::unique_ptr<itkMeshOperations > m_meshOperations;
        std::unique_ptr<sampling::MarkovChain<mhfitting::MHFittingParameters> > m_chain; // FIXME change type
        PreprocessedImagePointerType m_preprocessedTargetImage;
        SamplerPointerType m_sampler;

        std::unique_ptr<::sampling::BestMatchLogger<mhfitting::MHFittingParameters> > m_bestMatchLogger;
        std::unique_ptr<::sampling::LoggerDistributor<mhfitting::MHFittingParameters> > m_loggerDistributor;
        std::unique_ptr<::sampling::AcceptRejectLogger<mhfitting::MHFittingParameters> > m_printLogger;
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
