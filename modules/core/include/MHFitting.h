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
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANYlineGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#ifndef STATISMO_MHFITTING_H
#define STATISMO_MHFITTING_H

#include "ActiveShapeModel.h"
#include "ASMFitting.h"
#include "ASMPointSampler.h"
#include "MultiVariateNormalDistribution.h"
#include <boost/thread.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread/future.hpp>

#include <Sampling/DistributionEvaluator.h>
#include <Sampling/Evaluator/ProductEvaluator.h>
#include <Sampling/MarkovChain.h>
#include <Sampling/Algorithm/ForwardChain.h>
#include <Sampling/Algorithm/Metropolis.h>
#include <Sampling/RandomGenerator.h>
#include <Sampling/ProposalGenerator.h>
#include <Sampling/Proposal/RandomProposal.h>
#include <Sampling/Logger/QuietLogger.h>
#include <Sampling/RandomGenerator.h>
#include <Sampling/Algorithm/MetropolisHastings.h>
#include <Sampling/Logger/BestMatchLogger.h>
#include <Sampling/Proposal/MarkovChainProposal.h>


namespace statismo {
  using namespace sampling;
  using namespace std;





    class MHFittingConfiguration {
    private:

    public:
        MHFittingConfiguration(const ASMFittingConfiguration& asmFittingConfiguration, unsigned numberOfProfilePoints = 100 ) : m_asmFittingConfiguration(asmFittingConfiguration) { }
        const ASMFittingConfiguration& GetAsmFittingconfiguration() const {return m_asmFittingConfiguration; }

    private:
        ASMFittingConfiguration m_asmFittingConfiguration;
    };







  class MHFittingParameters {
    public:
      MHFittingParameters() {
      }

      MHFittingParameters(VectorType coefficients, VectorType rigidParameters) :
        m_coefficients(coefficients),
        m_rigidParameters(rigidParameters) { }

      VectorType GetCoefficients() const {
        return m_coefficients;
      }

      VectorType GetRigidTransformParameters() const {
        return m_rigidParameters;
      }


      int size() {
        // TODO: if this is really used we should change to a different vectorized representation of all parameters
        return m_coefficients.size() + m_rigidParameters.size();
      }

      float& operator[]( int i ) {
        // TODO: if this is really used we should change to a different vectorized representation of all parameters
        throw new std::runtime_error("operator[] not implemented for MHFittingResult");
        return m_coefficients(0,0);
      }

    private:
      VectorType m_coefficients;
        VectorType m_rigidParameters;
  };



    template <class MeshType, class PointType>
    class MeshOperations {
    public:
        virtual std::pair<PointType, long> findClosestPoint(MeshType mesh, PointType pt) const = 0;
        virtual VectorType normalAtPoint(MeshType mesh, PointType pt) const = 0;
        virtual short huAtPoint(MeshType mesh, long ptId) const = 0;
        virtual MeshType transformMesh(const MHFittingParameters& fittingParameters) const = 0;
        virtual PointType getPointWithId(MeshType mesh, unsigned id) const = 0;
        virtual PointType transformToModelSpace(const VectorType& rigidTransformParameters, PointType pt) const = 0;
    };



  template< class TPointSet, class TImage>
  class mcmc {

  public:

      typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
      typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
      typedef typename Representer<TPointSet>::PointType PointType;
      typedef StatisticalModel<TPointSet> StatisticalModelType;
      typedef typename StatisticalModelType::RepresenterType RepresenterType;
      typedef typename StatisticalModelType::DatasetPointerType SampleType;
      typedef typename StatisticalModelType::ValueType ValueType;
      typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
      typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
      typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
      typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;
      typedef std::vector<std::pair<unsigned, PointType> > CorrespondencePoints;



      template <class T>
      class UncertaintyLogger : public ChainLogger<MHFittingParameters> {
          typedef typename Representer<T>::PointType PointType;
          typedef unsigned int PointId;

      public:
          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> MeshOperationsType;
          UncertaintyLogger(const Representer<T>* representer, const MeshOperationsType* meshOps, const CorrespondencePoints& correspondencePoints,  const vector< PointType >& targetPoints)

                  : m_representer(representer),
                    m_meshOps(meshOps),
                    m_correspondencePoints(correspondencePoints),
                    m_targetPoints(targetPoints)
          {
              // fill the map with an empty vector of points
              for( unsigned i = 0; i < m_correspondencePoints.size(); ++i)  {
                  unsigned id = m_correspondencePoints[i].first;
                  m_uncertaintyPerCorrespondingPoint.insert(std::make_pair(id, std::vector<PointType>()));
              }


          }
          ~UncertaintyLogger() {
              std::cout << "destructor of uncertainty logger" << std::endl;

          }


          virtual void notifyAccept(
                  const MHFittingParameters& parameters,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {

              typename Representer<T>::DatasetPointerType sample = m_meshOps->transformMesh(parameters);
              for( unsigned i = 0; i < m_correspondencePoints.size(); ++i)  {
                  unsigned id = m_correspondencePoints[i].first;
                  PointType pointOnSample = m_meshOps->getPointWithId(sample, id);
                  m_uncertaintyPerCorrespondingPoint[id].push_back(pointOnSample);
              }

              // we keep the parameters, such that we can compute the average at the end
              m_params.push_back(parameters);
          }

          /** \brief Function gets called whenever an algorithm rejects a proposal */
          virtual void notifyReject(
                  const MHFittingParameters& sample,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {}

          /** \brief Function gets called whenever an algorithm is reset to a new state */
          virtual void notifyReset (
                  const MHFittingParameters& state,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {}


          MultiVariateNormalDistribution uncertaintyForCorrespondencePoint(unsigned id) {

              if (m_params.size() == 0) {
                  throw std::runtime_error("we did not have an accepted sample in the unceratinty estimation");
              }
              VectorType avgTransformsParameters = VectorType::Zero(m_params.front().GetRigidTransformParameters().rows());
              for (typename std::list<MHFittingParameters>::const_iterator it = m_params.begin(); it != m_params.end(); ++it) {
                  avgTransformsParameters += it->GetRigidTransformParameters();
              }
              avgTransformsParameters /= m_params.size();

              VectorType mean = VectorType::Zero(3);
              MatrixType cov = MatrixType::Zero(3, 3);


              const std::vector<PointType>& v = m_uncertaintyPerCorrespondingPoint[id];


              for (typename std::vector<PointType>::const_iterator it = v.begin(); it != v.end(); ++it) {
                  VectorType ptAsVec =  m_representer->PointToVector(m_meshOps->transformToModelSpace(avgTransformsParameters, *it));
                  mean += ptAsVec;
              }
              mean /= v.size();


              for (typename std::vector<PointType>::const_iterator it = v.begin(); it != v.end(); ++it) {
                  VectorType ptAsVec =  m_representer->PointToVector(m_meshOps->transformToModelSpace(avgTransformsParameters, *it));
                  cov += (ptAsVec - mean) * (ptAsVec -mean).transpose();
              }
              cov /= (v.size() - 1 );

              MultiVariateNormalDistribution mvn(mean, cov);
              return mvn;

          }

      private:
          const Representer<T>* m_representer;
          const MeshOperationsType* m_meshOps;
          std::vector<PointType> m_targetPoints;
          CorrespondencePoints m_correspondencePoints;
          std::list<MHFittingParameters> m_params;
          std::map<PointId, std::vector<PointType> > m_uncertaintyPerCorrespondingPoint;
      };



      template <class T>
      class InLungLogger : public ChainLogger<MHFittingParameters> {
          typedef typename Representer<T>::PointType PointType;
          typedef unsigned int PointId;

      public:


          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> MeshOperationsType;
          InLungLogger(ChainLogger<MHFittingParameters>* otherLogger, const Representer<T>* representer, const MeshOperationsType* meshOps, const char* loggerName)

                  : m_compositeLogger(otherLogger),
                    m_representer(representer),
                    m_meshOps(meshOps),
                    m_loggerName(loggerName)
          {
          }
          ~InLungLogger() {
              std::cout << "destructor of uncertainty logger" << std::endl;

          }


          virtual void notifyAccept(
                  const MHFittingParameters& parameters,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {


              m_compositeLogger->notifyAccept(parameters, dProbValue, proposal, evaluator);


              typename Representer<T>::DatasetPointerType sample = m_meshOps->transformMesh(parameters);
              unsigned numPoints = m_representer->GetDomain().GetNumberOfPoints();
              bool isInsideBone = false;
              unsigned numInBone = 0;
              for( int i = 0; i <numPoints; ++i) {

                  if (m_meshOps->huAtPoint(sample, i) > 1150 || m_meshOps->huAtPoint(sample, i) < 400 )   {
//                      isInsideBone = true;
//                      break;
                      numInBone +=  1;
                  }
              }

//              if (isInsideBone) {
//                  distance = -std::numeric_limits<double>::infinity();
//
//              } else {
//                  distance = 0;
//              }

              std::cout << "Accepted Step in " << m_loggerName << " : New number in bone is " << numInBone << " parameter norm " << parameters.GetCoefficients().squaredNorm() << std::endl;

          }

          /** \brief Function gets called whenever an algorithm rejects a proposal */
          virtual void notifyReject(
                  const MHFittingParameters& sample,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {
              m_compositeLogger->notifyReject(sample, dProbValue, proposal, evaluator);

          }

          /** \brief Function gets called whenever an algorithm is reset to a new state */
          virtual void notifyReset (
                  const MHFittingParameters& state,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {
              m_compositeLogger->notifyReset(state, dProbValue, proposal, evaluator);
          }




      private:
          ChainLogger<MHFittingParameters>* m_compositeLogger;
          const Representer<T>* m_representer;
          const MeshOperationsType* m_meshOps;
          const char* m_loggerName;
      };



      template <class T>
      class PoseLogger : public ChainLogger<MHFittingParameters> {

      public:
          PoseLogger(const char* loggerName)
:                    m_loggerName(loggerName)
          {
          }


          virtual void notifyAccept(
                  const MHFittingParameters& parameters,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {



              std::cout << "accepted propseal in logger " << m_loggerName << "with parameters " << parameters.GetRigidTransformParameters().transpose()  << std::endl;


          }

          /** \brief Function gets called whenever an algorithm rejects a proposal */
          virtual void notifyReject(
                  const MHFittingParameters& sample,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {

              std::cout << "rejected propseal in logger " << m_loggerName << "with parameters " << sample.GetRigidTransformParameters().transpose() << std::endl;
          }

          /** \brief Function gets called whenever an algorithm is reset to a new state */
          virtual void notifyReset (
                  const MHFittingParameters& state,
                  const double& dProbValue,
                  ProposalGeneratorInterface<MHFittingParameters>* proposal,
                  DistributionEvaluatorInterface<MHFittingParameters>* evaluator
          ) {}



      private:
          const char* m_loggerName;
      };




      // Proposals
      /*============================================================================*/
//      class GaussianModelUpdate : public ProposalGenerator<MHFittingParameters > {
//        private:
//          double sigmaShape;
//          double sigmaTranslation;
//          double sigmaRotation;
//          unsigned m_maxNumberOfShapeParameters;
//          RandomGenerator* rgen;
//
//        public:
//          GaussianModelUpdate( double stepSizeShape, double stepSizeRotation, double stepSizeTranslation, unsigned maxNumberOfShapeParameters, RandomGenerator* rgen )
//                  : sigmaShape(stepSizeShape), sigmaRotation(stepSizeRotation), sigmaTranslation(stepSizeTranslation), rgen(rgen), m_maxNumberOfShapeParameters(maxNumberOfShapeParameters) {}
//
//          // ProposalGeneratorInterface interface
//        public:
//          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){
//            VectorType shapeParams = currentSample.GetCoefficients();
//
//              VectorType newShapeParams(shapeParams.size());
//              newShapeParams.fill(0);
//            for (unsigned i = 0; i < std::min(m_maxNumberOfShapeParameters, static_cast<unsigned>(shapeParams.size())); ++i) {
//                newShapeParams[i] = shapeParams[i] + rgen->normalDbl() * sigmaShape;
//            }
//
//
//              VectorType newRigidParams = currentSample.GetRigidTransformParameters();
//              for (unsigned i = 0; i < currentSample.GetRigidTransformParameters().size(); ++i) {
//                  if (i < 3) {// rotation parameters
//                    newRigidParams[i] += rgen->normalDbl() * sigmaRotation * 0.01;
//                  }
//                  else { // tranlation scale
//                      newRigidParams[i] += rgen->normalDbl() * sigmaTranslation;
//                  }
//              }
//            proposal = MHFittingParameters(newShapeParams, newRigidParams);
//          }
//
//          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
//            return 0.0;
//          }
//      };
      class RandomShapeUpdate : public ProposalGenerator<MHFittingParameters > {
        private:
          double sigmaShape;
          double sigmaTranslation;
          double sigmaRotation;
          unsigned m_maxNumberOfShapeParameters;
          RandomGenerator* rgen;

        public:
          RandomShapeUpdate( double stepSizeShape, unsigned maxNumberOfShapeParameters, RandomGenerator* rgen )
                  : sigmaShape(stepSizeShape),  rgen(rgen), m_maxNumberOfShapeParameters(maxNumberOfShapeParameters) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){
            VectorType shapeParams = currentSample.GetCoefficients();

              VectorType newShapeParams(shapeParams.size());
              newShapeParams.fill(0);
            for (unsigned i = 0; i < std::min(m_maxNumberOfShapeParameters, static_cast<unsigned>(shapeParams.size())); ++i) {
                newShapeParams[i] = shapeParams[i] + rgen->normalDbl() * sigmaShape;
            }


            proposal = MHFittingParameters(newShapeParams, currentSample.GetRigidTransformParameters());
          }

          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
          }
      };


        class RotationUpdate : public ProposalGenerator<MHFittingParameters > {
        private:
            unsigned m_axis;
          double m_sigmaRotation;
          RandomGenerator* m_rgen;

        public:
            RotationUpdate(unsigned axis, double stepSizeRotation, RandomGenerator* rgen)
                  : m_axis(axis) , m_sigmaRotation(stepSizeRotation),  m_rgen(rgen) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

              VectorType newRigidParams = currentSample.GetRigidTransformParameters();
              newRigidParams[m_axis] += m_rgen->normalDbl() * m_sigmaRotation;
            proposal = MHFittingParameters( currentSample.GetCoefficients(), newRigidParams);
          }

          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
          }
      };


      class TranslationUpdate : public ProposalGenerator<MHFittingParameters > {
      private:
          unsigned m_axis;
          double m_sigmaTranslation;
          RandomGenerator* m_rgen;

      public:
          TranslationUpdate(unsigned axis, double stepSizeTranslations, RandomGenerator* rgen)
                  : m_axis(axis) , m_sigmaTranslation(stepSizeTranslations),  m_rgen(rgen) {}

          // ProposalGeneratorInterface interface
      public:
          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

              VectorType newRigidParams = currentSample.GetRigidTransformParameters();
              newRigidParams[3 + m_axis] += m_rgen->normalDbl() * m_sigmaTranslation;
              proposal = MHFittingParameters( currentSample.GetCoefficients(), newRigidParams);
          }

          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
              return 0.0;
          }
      };


      /*
      class ASMModelUpdate : public ProposalGenerator<MHFittingParameters > {
        private:
          int N;
          ActiveShapeModelType* m_asmodel;
          ASMFittingConfiguration m_fittingConfiguration;
          PreprocessedImageType* m_target;
          PointSamplerType* m_sampler;

        public:
          ASMModelUpdate(ASMFittingConfiguration config, ActiveShapeModelType* asmodel, PreprocessedImageType* target,  PointSamplerType* sampler, int nSteps = 1 )
                  : m_fittingConfiguration(config),
                  m_asmodel(asmodel) ,
                    m_target(target),
                    m_sampler(sampler),
                  N(nSteps) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_fittingConfiguration, m_asmodel, currentSample.GetCoefficients(), currentSample.GetRigidTransform(), m_target, m_sampler);
            ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();
              statismo::VectorType newCoeffs = result.GetCoefficients();
              statismo::VectorType currCoeffs = currentSample.GetCoefficients();
              statismo::VectorType newProposal = currCoeffs + (newCoeffs - currCoeffs) * 0.1;
            proposal = MHFittingParameters(newProposal,result.GetRigidTransform());
          }
          virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
          }
      };
       */

      // Evaluators
      /*============================================================================*/
      class PositionEvaluator : public DistributionEvaluator<PointType>
      {
        public:
          /// dummy gradient, zero, do nothing
          virtual void evalGradient( PointType& gradient, PointType const& sample )
          {
            throw logic_error("no gradient implemented in PositionEvaluator");
          }

          /// Returns whether this evaluator is separable into x and y coordinates
          virtual bool isSeparable() const
          {
            return false;
          }
      };


      class ModelPriorEvaluator : public DistributionEvaluator<MHFittingParameters>
      {
      public:
          virtual double evalSample(const MHFittingParameters& currentSample) {

                VectorType coefficients = currentSample.GetCoefficients();

              return -0.5 * coefficients.size() *log( 2*M_PI ) - 0.5 * coefficients.squaredNorm();

          }

          /// dummy gradient, zero, do nothing
          virtual void evalGradient( PointType& gradient, PointType const& sample )
          {
              throw logic_error("no gradient implemented in PositionEvaluator");
          }

          /// Returns whether this evaluator is separable into x and y coordinates
          virtual bool isSeparable() const
          {
              return false;
          }
      };


      template <class T>
      class InLungOrBoneEvaluator : public DistributionEvaluator< MHFittingParameters > {
          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> ClosestPointType;

      public:
          InLungOrBoneEvaluator( const Representer<T>* representer, const ClosestPointType* closestPoint, ActiveShapeModelType* asmodel) :
                  m_asmodel(asmodel),
                  m_meshOperations(closestPoint)

          {}

          // DistributionEvaluatorInterface interface
      public:
          virtual double evalSample(const MHFittingParameters& currentSample) {
              double distance = 0.0;

              // TODO to draw full sample only for landmarks is inefficient
//              typename RepresenterType::DatasetPointerType  sampleShape = m_asmodel->GetStatisticalModel()->DrawSample(currentSample.GetCoefficients());
//              typename RepresenterType::DatasetPointerType sample = m_asmodel->GetRepresenter()->TransformMesh(sampleShape, currentSample.GetRigidTransform());


              //std::cout << currentSample.GetRigidTransform()->GetParameters() << std::endl;

              typename RepresenterType::DatasetPointerType sample = m_meshOperations->transformMesh(currentSample);
              unsigned numPoints = m_asmodel->GetRepresenter()->GetDomain().GetNumberOfPoints();
              bool isInsideBone = false;
              double increase = -std::numeric_limits<double>::max() / numPoints;
              unsigned numInBone = 0;
              for( int i = 0; i <numPoints; ++i) {
                  //statismo::VectorType normal = m_closestPoint->normalAtPoint(sample, m_tra)
                  //PointType closestPtOnSample =  m_closestPoint->findClosestPoint(sample, m_targetPoints[i]).first;
                  //double d = m_eval->evalSample(closestPtOnSample-m_targetPoints[i]);


                  if (m_meshOperations->huAtPoint(sample, i) > 1150 || m_meshOperations->huAtPoint(sample, i) < 400 )   {
//                      isInsideBone = true;
//                      break;

                      distance += increase;
                      numInBone +=  1;
                  }
                  else {
                      distance += 0;
                  }
              }

//              if (isInsideBone) {
//                  distance = -std::numeric_limits<double>::infinity();
//
//              } else {
//                  distance = 0;
//              }

              return distance;

          }
      private:
          const vector< PointType > m_targetPoints;
          const ActiveShapeModelType* m_asmodel;
          PositionEvaluator* m_eval;
          const ClosestPointType* m_meshOperations;

      };




      class Gaussian3DPositionDifferenceEvaluator : public PositionEvaluator
      {
        private:
          const RepresenterType* rep;
        public:
          Gaussian3DPositionDifferenceEvaluator( const RepresenterType* representer, double sigma = 1.0) : rep(representer){
            m_sigma = sigma;
            m_dNormalizer = -1.5*log( 2*M_PI ) - 3.0 * log( sigma );
          }


          double evalSample( PointType const& diff )
          {
            double dVal = rep->PointToVector(diff).squaredNorm(); // TODO check if PointType has this function norm2
            return -dVal/(2.0*m_sigma*m_sigma) + m_dNormalizer;
          }

          void evalGradient( PointType& grad, PointType const& diff )
          {
            //grad = diff*(-1) / (sigma * sigma);
          }

          virtual bool isSeparable() const
          {
            return true;
          }

        private:
          double m_sigma;
          double m_dNormalizer;
      };


      template <class T>
      class PointEvaluator : public DistributionEvaluator< MHFittingParameters > {
          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> ClosestPointType;

        public:
          PointEvaluator( const Representer<T>* representer, const ClosestPointType* closestPoint, const CorrespondencePoints& correspondencePoints, ActiveShapeModelType* asmodel, PositionEvaluator* evaluator) :
            m_correspondencePoints(correspondencePoints),
            m_asmodel(asmodel),
            m_eval(evaluator),
            m_closestPoint(closestPoint)

          {}

          // DistributionEvaluatorInterface interface
        public:
          virtual double evalSample(const MHFittingParameters& currentSample) {
            double distance = 0.0;

            // TODO to draw full sample only for landmarks is inefficient

              typename RepresenterType::DatasetPointerType sample = m_closestPoint->transformMesh(currentSample);


              for( unsigned i = 0; i < m_correspondencePoints.size(); ++i)  {

                  PointType pointOnSample = m_closestPoint->getPointWithId(sample, m_correspondencePoints[i].first);
                  double d = m_eval->evalSample(pointOnSample - m_correspondencePoints[i].second);
                  distance += d;
              }

            return distance;

          }
        private:
          CorrespondencePoints m_correspondencePoints;
          const ActiveShapeModelType* m_asmodel;
          PositionEvaluator* m_eval;
          const ClosestPointType* m_closestPoint;

      };

      template <class T>
      class LineEvaluator : public DistributionEvaluator< MHFittingParameters > {
          typedef MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> ClosestPointType;
          typedef std::map<float, std::vector<VectorType> > LineMapType;
      public:
          LineEvaluator( const Representer<T>* representer, const ClosestPointType* closestPoint, const vector< PointType >& targetPoints, ActiveShapeModelType* asmodel, double noisevar) :
                  m_targetPoints(targetPoints),
                  m_asmodel(asmodel),
                  m_closestPoint(closestPoint),
                  m_noiseVar(noisevar)
          {

              // TODO HACK, we separate the points into lines
              for (unsigned i = 0 ; i < targetPoints.size(); ++i ) {
                  VectorType pt = representer->PointToVector(targetPoints[i]);
                  if (m_lineMap.find(pt(2)) != m_lineMap.end()) {
                        m_lineMap.insert(std::make_pair(pt(2), std::vector<VectorType>()));
                  }
                  m_lineMap[pt(2)].push_back(pt);

              }
              std::cout << "number of lines " << m_lineMap.size();
              VectorType mean = VectorType::Zero(1);
              mean(0) = 0.0;
              MatrixType cov = MatrixType::Zero(1,1);
              cov(0,0)=m_noiseVar;
                m_likelihoodModel = MultiVariateNormalDistribution(mean, cov);
          }

          // DistributionEvaluatorInterface interface
      public:


          double likelihoodForLine(const std::vector<VectorType> & line, typename RepresenterType::DatasetPointerType sample) {
//              double sumOfsquaraedDistance = 0.0;
//              for( int i = 0; i < m_targetPoints.size(); ++i) {
//
//                  PointType closestPtOnSample =  m_closestPoint->findClosestPoint(sample, m_targetPoints[i]).first;
//                  double d = (closestPtOnSample-m_targetPoints[i]).GetNorm();
//
//                  sumOfsquaraedDistance += d * d;
//              }
//              double avgDistance = sumOfsquaraedDistance / m_targetPoints.size();
//              VectorType avgDistanceAsVec(1);
//              avgDistanceAsVec << avgDistance;
//              m_likelihoodModel.logpdf(avgDistanceAsVec);

              double sumOfLog = 0.0;
              for( int i = 0; i < m_targetPoints.size(); ++i) {

                  PointType closestPtOnSample =  m_closestPoint->findClosestPoint(sample, m_targetPoints[i]).first;
                  double d = (closestPtOnSample-m_targetPoints[i]).GetNorm();
                  VectorType distAsVec(1);
                  distAsVec << d;
                  sumOfLog += m_likelihoodModel.logpdf(distAsVec);
              }
                return sumOfLog;


          }

          virtual double evalSample(const MHFittingParameters& currentSample) {


              // TODO to draw full sample only for landmarks is inefficient

              typename RepresenterType::DatasetPointerType sample = m_closestPoint->transformMesh(currentSample);

              double sumLikelihood = 0;
              for (LineMapType::iterator it = m_lineMap.begin(); it != m_lineMap.end(); ++it) {
                  sumLikelihood += likelihoodForLine(it->second, sample);
              }
             return sumLikelihood;
          }
      private:
          MultiVariateNormalDistribution m_likelihoodModel;
          const vector< PointType > m_targetPoints;
          const ActiveShapeModelType* m_asmodel;
          const ClosestPointType* m_closestPoint;
          LineMapType m_lineMap;
          double m_noiseVar;
      };



      struct ASMLikelihoodForChunk {


          ASMLikelihoodForChunk(double _aggregatedLikelihood): aggregatedLikelihood(_aggregatedLikelihood) {  }

          double aggregatedLikelihood;

          // emulate move semantics, as boost::async seems to depend on it.
          ASMLikelihoodForChunk& operator=(BOOST_COPY_ASSIGN_REF(ASMLikelihoodForChunk) rhs) { // Copy assignment

              if (&rhs != this) {
                  copyMembers(rhs);
              }
              return *this;
          }

          ASMLikelihoodForChunk(BOOST_RV_REF(ASMLikelihoodForChunk) that) { //Move constructor
              copyMembers(that);
          }
          ASMLikelihoodForChunk& operator=(BOOST_RV_REF(ASMLikelihoodForChunk) rhs) { //Move assignment
              if (&rhs != this) {
                  copyMembers(rhs);
              }
              return *this;
          }
      private:
      BOOST_COPYABLE_AND_MOVABLE(ASMLikelihoodForChunk)
          void copyMembers(const ASMLikelihoodForChunk& that) {
              aggregatedLikelihood = that.aggregatedLikelihood;
          }
      };



/*
      template <class T>
      class ASMEvaluator : public DistributionEvaluator< MHFittingParameters > {
        typedef  ASMPreprocessedImage<typename RepresenterType::DatasetType> PreprocessedImageType;

      public:
          ASMEvaluator(ActiveShapeModelType* asmodel, PreprocessedImageType* image, PointSamplerType* sampler) :
                  m_asmodel(asmodel),
                  m_image(image),
                  m_sampler(sampler)
          {}

          // DistributionEvaluatorInterface interface
      public:
          virtual double evalSample(const MHFittingParameters& currentSample) {


              typename RepresenterType::DatasetPointerType sample = m_closestPoint->TransformMesh(m_asmodel->GetStatisticalModel(), currentSample);


              unsigned numProfilePointsUsed = 500;
              unsigned step = m_asmodel->GetProfiles().size() / numProfilePointsUsed;
//              FeatureExtractorType* fe = m_asmodel->GetFeatureExtractor()->CloneForTarget(m_asmodel,currentSample.GetCoefficients(),currentSample.GetRigidTransform());

              unsigned numChunks =  boost::thread::hardware_concurrency() + 1;
              std::vector<boost::future<ASMLikelihoodForChunk>* > futvec;


              for (unsigned i = 0; i < numChunks; ++i) {
                  unsigned profileIdStart = m_asmodel->GetProfiles().size() / numChunks  * i;
                  unsigned profileIdEnd = m_asmodel->GetProfiles().size() / numChunks * (i + 1);

                  boost::future<ASMLikelihoodForChunk> *fut = new boost::future<ASMLikelihoodForChunk>(
                          boost::async(boost::launch::async, boost::bind(&ASMEvaluator<T>::evalSampleForProfiles,
                                                                         this, profileIdStart, profileIdEnd, step,// fe,
                                                                         currentModelInstance, currentSample)));
                  futvec.push_back(fut);
              }


              double loglikelihood = 0.0;
              for (unsigned i = 0; i < futvec.size(); i++) {

                  ASMLikelihoodForChunk likelihoodForChunk = futvec[i]->get();
                  loglikelihood += likelihoodForChunk.aggregatedLikelihood;
                  delete futvec[i];
              }


              return loglikelihood;

          }


          ASMLikelihoodForChunk evalSampleForProfiles(unsigned profileIdStart, unsigned profileIdEnd, unsigned step, typename RepresenterType::DatasetPointerType currentModelInstance, const MHFittingParameters& currentSample) {

              FeatureExtractorType* fe = m_asmodel->GetFeatureExtractor()->CloneForTarget(m_asmodel,currentSample.GetCoefficients(),currentSample.GetRigidTransform());

              double loglikelihood = 0.0;

              for (unsigned i = profileIdStart; i < profileIdEnd; i += step) {

                  ASMProfile profile = m_asmodel->GetProfiles()[i];
                  long ptId = profile.GetPointId();

                  statismo::VectorType features;


                  bool ok = fe->ExtractFeatures(features, m_image, currentModelInstance->GetPoint(ptId));

                  if (ok) {
                      loglikelihood += profile.GetDistribution().logpdf(features);
                  } else {
                      std::cout << "feasutre not ok " << std::endl;
                  }
              }
              // evaluate profile points at ...
              fe->Delete();
//              ASMLikelihoodForChunk asmLikelihoodForChunk(loglikelihood);
              return ASMLikelihoodForChunk(loglikelihood);
          }

      private:
          const ActiveShapeModelType* m_asmodel;
          PreprocessedImageType* m_image;

          PointSamplerType* m_sampler;
          RandomGenerator* m_rGen;
      };
*/

      // TODO: Full image evaluator, i.e. ASM Evaluator is missing


    public:
      // "Script"
      /*============================================================================*/
      template <class T>
      class BasicSampling {
        public:






          static MarkovChain<MHFittingParameters >* buildLmAndHuChain(
            const Representer<T>* representer,
            const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* meshOperations,
            CorrespondencePoints correspondencePoints,
            vector<PointType> targetPoints,
            ActiveShapeModelType* asmodel,
            MHFittingParameters& initialParameters,
            ChainLogger<MHFittingParameters>* logger ) {

              unsigned numPCAComponents = asmodel->GetStatisticalModel()->GetNumberOfPrincipalComponents();

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);
            MHFittingParameters init = MHFittingParameters(initialParameters.GetCoefficients(), initialParameters.GetRigidTransformParameters());

              RandomShapeUpdate* shapeUpdateRough = new RandomShapeUpdate(0.2, numPCAComponents, rGen);
              RandomShapeUpdate* shapeUpdateFine = new RandomShapeUpdate(0.1, numPCAComponents, rGen);
              RandomShapeUpdate* shapeUpdateFinest = new RandomShapeUpdate(0.025, numPCAComponents, rGen);
              RotationUpdate* rotUpdateX = new RotationUpdate(0, 0.01, rGen);
              RotationUpdate* rotUpdateY = new RotationUpdate(1, 0.01, rGen);
              RotationUpdate* rotUpdateZ = new RotationUpdate(2, 0.01, rGen);

              vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> rotUpdateVec(3);
              rotUpdateVec[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateX, 0.8);
              rotUpdateVec[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateY, 0.1);
              rotUpdateVec[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateZ, 0.1);
              RandomProposal<MHFittingParameters >* rotUpdate = new RandomProposal<MHFittingParameters >(rotUpdateVec, rGen);

              // TODO Note: there is currently no random translation

            vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector(4);
            gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateRough, 0.1);
            gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateFine, 0.2);
              gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateFinest, 0.4);
            gaussMixtureProposalVector[3] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdate, 0.2);


            RandomProposal<MHFittingParameters >* gaussMixtureProposal = new RandomProposal<MHFittingParameters >(gaussMixtureProposalVector, rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 1.0);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, meshOperations, correspondencePoints, asmodel, diffEval);
              LineEvaluator<T>* lineEval = new LineEvaluator<T>(representer, meshOperations, targetPoints, asmodel, 0.1);

            ModelPriorEvaluator* modelPriorEvaluator = new ModelPriorEvaluator();

            InLungOrBoneEvaluator<T>* huEvaluator = new InLungOrBoneEvaluator<T>(representer, meshOperations, asmodel);

            std::vector<DistributionEvaluator<MHFittingParameters >*> huEvaluatorList;
              //  huEvaluatorList.push_back(huEvaluator);
            //huEvaluatorList.push_back(modelPriorEvaluator);
//              huEvaluatorList.push_back(pointEval);
              huEvaluatorList.push_back(lineEval);


              QuietLogger<MHFittingParameters>* ql = new QuietLogger<MHFittingParameters>();
              InLungLogger <T>* loggerFilterChain = new InLungLogger<T>(logger, representer, meshOperations, "filter chain");
            MarkovChain<MHFittingParameters >* huChain = new MetropolisHastings<MHFittingParameters >(gaussMixtureProposal, new ProductEvaluator<MHFittingParameters>(huEvaluatorList), ql, init, rGen);
            MarkovChainProposal<MHFittingParameters>* huChainProposal = new MarkovChainProposal<MHFittingParameters>(huChain, 1, true);


            std::vector<DistributionEvaluator<MHFittingParameters >*> lmAndHuEvaluatorList;
                //lmAndHuEvaluatorList.push_back(pointEval);
              //lmAndHuEvaluatorList.push_back(lineEval);
              lmAndHuEvaluatorList.push_back(huEvaluator);
              lmAndHuEvaluatorList.push_back(modelPriorEvaluator);

              //InLungLogger <T>* loggerFinalChain = new InLungLogger<T>(representer, meshOperations, "final chain");
            MarkovChain<MHFittingParameters >* lmAndHuChain = new MetropolisHastings<MHFittingParameters >(huChainProposal, new ProductEvaluator<MHFittingParameters>(lmAndHuEvaluatorList), loggerFilterChain, init, rGen);
            return lmAndHuChain;
          }




          // estimate the uncertainty at the given correspondence Points, by sampling from the initialPoseChain.
          static std::map<unsigned, MultiVariateNormalDistribution> estimatePointUncertaintyForInitialPoseChain(const Representer <T> *representer,
                                                                  const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> *meshOperations,
                                                                  CorrespondencePoints correspondencePoints,
                                                                  vector<PointType> targetPoints,
                                                                  ActiveShapeModelType *asmodel,
                                                                  MHFittingParameters &initialParameters)
          {
              UncertaintyLogger<T>* ul = new UncertaintyLogger<T>(representer, meshOperations, correspondencePoints, targetPoints);
              MarkovChain<MHFittingParameters>* chain = buildInitialPoseChain(representer, meshOperations, correspondencePoints, targetPoints, asmodel, initialParameters, ul);

              MHFittingParameters params;//(m_sourceTransform,m_sourceCoefficients);
              for (unsigned i = 0; i < 1000; ++i) {
                  chain->next(params);
              }

              std::map<unsigned, MultiVariateNormalDistribution> uncertaintyMap;

              for (unsigned i = 0 ; i < correspondencePoints.size(); ++i) {
                  unsigned id = correspondencePoints[i].first;

                  uncertaintyMap.insert(std::make_pair(id, ul->uncertaintyForCorrespondencePoint(id)));
              }
              return uncertaintyMap;

          }

          static MarkovChain<MHFittingParameters >* buildInitialPoseChain(
            const Representer<T>* representer,
            const MeshOperations<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* closestPoint,
            CorrespondencePoints correspondencePoints,
            vector<PointType> targetPoints,
            ActiveShapeModelType* asmodel,
            MHFittingParameters& initialParameters,
            ChainLogger<MHFittingParameters>* logger ) {

              unsigned numPCAComponents = asmodel->GetStatisticalModel()->GetNumberOfPrincipalComponents();

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);


              RotationUpdate* rotUpdateX = new RotationUpdate(0, 0.01, rGen);
              RotationUpdate* rotUpdateY = new RotationUpdate(1, 0.01, rGen);
              RotationUpdate* rotUpdateZ = new RotationUpdate(2, 0.01, rGen);

              vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> rotUpdateVec(3);
              rotUpdateVec[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateX, 0.6);
              rotUpdateVec[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateY, 0.2);
              rotUpdateVec[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateZ, 0.2);
              RandomProposal<MHFittingParameters >* rotUpdate = new RandomProposal<MHFittingParameters >(rotUpdateVec, rGen);


              TranslationUpdate* transUpdateXSmall = new TranslationUpdate(0, 1, rGen);
              TranslationUpdate* transUpdateYSmall = new TranslationUpdate(1, 1, rGen);
              TranslationUpdate* transUpdateZSmall = new TranslationUpdate(2, 1, rGen);
              TranslationUpdate* transUpdateXLarge = new TranslationUpdate(0, 10, rGen);
              TranslationUpdate* transUpdateYLarge = new TranslationUpdate(1, 10, rGen);
              TranslationUpdate* transUpdateZLarge = new TranslationUpdate(2, 10, rGen);


              vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> transUpdateVec(6);
              transUpdateVec[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateXSmall, 0.1);
              transUpdateVec[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateYSmall, 0.1);
              transUpdateVec[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateZSmall, 0.1);
              transUpdateVec[3] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateXLarge, 0.1);
              transUpdateVec[4] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateYLarge, 0.1);
              transUpdateVec[5] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateZLarge, 0.1);

              RandomProposal<MHFittingParameters >* transUpdate = new RandomProposal<MHFittingParameters >(transUpdateVec, rGen);


              // TODO Note: there is currently no random translation
              RandomShapeUpdate* shapeUpdateRough = new RandomShapeUpdate(0.1, numPCAComponents, rGen);
              RandomShapeUpdate* shapeUpdateFine = new RandomShapeUpdate(0.05, numPCAComponents, rGen);
              RandomShapeUpdate* shapeUpdateFinest = new RandomShapeUpdate(0.01, numPCAComponents, rGen);

              vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector(5);

              gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateRough, 0.1);
              gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateFine, 0.1);
              gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateFinest, 0.1);
              gaussMixtureProposalVector[3] = pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdate, 0.1);
              gaussMixtureProposalVector[4] = pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdate, 0.1);

            RandomProposal<MHFittingParameters >* gaussMixtureProposal = new RandomProposal<MHFittingParameters >(gaussMixtureProposalVector, rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 1.0);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, closestPoint, correspondencePoints, asmodel, diffEval);
              LineEvaluator<T>* lineEval = new LineEvaluator<T>(representer, closestPoint, targetPoints, asmodel, 3.0);

              ModelPriorEvaluator* modelPriorEvaluator = new ModelPriorEvaluator();




            std::vector<DistributionEvaluator<MHFittingParameters >*> evaluatorList;
             evaluatorList.push_back(pointEval);
              evaluatorList.push_back(lineEval);
             evaluatorList.push_back(modelPriorEvaluator);

            MarkovChain<MHFittingParameters >* lmChain = new MetropolisHastings<MHFittingParameters >(gaussMixtureProposal, new ProductEvaluator<MHFittingParameters>(evaluatorList), logger, initialParameters, rGen);

            return lmChain;
          }


      };

  };





  template<typename TPointSet, typename TImage>
  class MHFittingStep {
      typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
      typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
      typedef typename Representer<TPointSet>::PointType PointType;
      typedef StatisticalModel<TPointSet> StatisticalModelType;
      typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
      typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
      typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
      typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;

      MHFittingStep(MarkovChain<MHFittingParameters >* chain)
        :
          m_chain(chain) { }

    private:
      class ProfileResult {
        public:
          unsigned int pointId;
          PointType candidatePoint;
          PointType transformedCandidatePoint;
      };


    public:


      static MHFittingStep *Create(MarkovChain<MHFittingParameters >* chain) {
        return new MHFittingStep( chain);
      }

      ~MHFittingStep() {
      }

      MHFittingParameters Perform() const {


        //          ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
        //          ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();


        // runs a markov chain and returns only the accepted proposals
        MHFittingParameters params;//(m_sourceTransform,m_sourceCoefficients);
        m_chain->next(params);
        MHFittingParameters mhResult(params.GetCoefficients(), params.GetRigidTransformParameters());
        return mhResult;
      }

    private:
      MarkovChain<MHFittingParameters >* m_chain;




  };
}

#endif //STATISMO_ASMFITTING_H


