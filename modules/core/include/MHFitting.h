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


#ifndef STATISMO_MHFITTING_H
#define STATISMO_MHFITTING_H

#include "ActiveShapeModel.h"
#include "ASMFitting.h"
#include "ASMPointSampler.h"
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




    template <class MeshType, class PointType>
    class ClosestPoint {
    public:
        virtual PointType findClosestPoint(MeshType mesh, PointType pt) const = 0;
    };

  template<typename RigidTransformPointerType>
  class MHFittingResult {
    public:
      MHFittingResult() {
      }

      MHFittingResult(VectorType coefficients, RigidTransformPointerType rigidTransform) :
        m_coefficients(coefficients),
        m_rigidTransform(rigidTransform) { }

      VectorType GetCoefficients() const {
        return m_coefficients;
      }

      RigidTransformPointerType GetRigidTransform() const {
        return m_rigidTransform;
      }

      int size() {
        // TODO: if this is really used we should change to a different vectorized representation of all parameters
        return m_coefficients.rows() * m_coefficients.cols() + m_rigidTransform->ParametersDimension;
      }

      float& operator[]( int i ) {
        // TODO: if this is really used we should change to a different vectorized representation of all parameters
        throw new std::runtime_error("operator[] not implemented for MHFittingResult");
        return m_coefficients(0,0);
      }

    private:
      VectorType m_coefficients;
      RigidTransformPointerType m_rigidTransform;
  };


  template< class TPointSet, class TImage>
  class mcmc {


      typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
      typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
      typedef TImage ImageType;
      typedef typename Representer<TPointSet>::PointType PointType;
      typedef StatisticalModel<TPointSet> StatisticalModelType;
      typedef typename StatisticalModelType::RepresenterType RepresenterType;
      typedef typename StatisticalModelType::DatasetPointerType SampleType;
      typedef typename StatisticalModelType::ValueType ValueType;
      typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
      typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
      typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
      typedef  ASMFittingStep<TPointSet, TImage> ASMFittingStepType;
      typedef MHFittingResult<RigidTransformPointerType> ChainSampleType;


      // Proposals
      /*============================================================================*/
      class GaussianModelUpdate : public ProposalGenerator<ChainSampleType > {
        private:
          double sigma;
          RandomGenerator* rgen;

        public:
          GaussianModelUpdate( double stepSize, RandomGenerator* rgen ) : sigma(stepSize), rgen(rgen) {}

          // ProposalGeneratorInterface interface
        public:
          virtual void generateProposal(ChainSampleType& proposal, const ChainSampleType& currentSample){
            VectorType shapeParams = currentSample.GetCoefficients();

              VectorType newParams(shapeParams.size());
            for (unsigned i = 0; i < shapeParams.size(); ++i) {
                newParams[i] = shapeParams[i] + rgen->normalDbl() * sigma;
            }

            proposal = ChainSampleType(newParams,currentSample.GetRigidTransform());
          }
          virtual double transitionProbability(const ChainSampleType& start, const ChainSampleType& end){
            return 0.0;
          }
      };


      class ASMModelUpdate : public ProposalGenerator<ChainSampleType > {
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
          virtual void generateProposal(ChainSampleType& proposal, const ChainSampleType& currentSample){

            ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_fittingConfiguration, m_asmodel, currentSample.GetCoefficients(), currentSample.GetRigidTransform(), m_target, m_sampler);
            ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();
	    
            proposal = ChainSampleType(result.GetCoefficients(),result.GetRigidTransform());
          }
          virtual double transitionProbability(const ChainSampleType& start, const ChainSampleType& end){
            return 0.0;
          }
      };


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


      class ModelPriorEvaluator : public DistributionEvaluator<ChainSampleType>
      {
      public:
          virtual double evalSample(const ChainSampleType& currentSample) {

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
      class PointEvaluator : public DistributionEvaluator< ChainSampleType > {
          typedef ClosestPoint<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType> ClosestPointType;

        public:
          PointEvaluator( const Representer<T>* representer, const ClosestPointType* closestPoint, const vector< PointType >& targetPoints, ActiveShapeModelType* asmodel, PositionEvaluator* evaluator) :
            m_targetPoints(targetPoints),
            m_asmodel(asmodel),
            m_eval(evaluator),
            m_closestPoint(closestPoint)

          {}

          // DistributionEvaluatorInterface interface
        public:
          virtual double evalSample(const ChainSampleType& currentSample) {
            double distance = 0.0;

            // TODO to draw full sample only for landmarks is inefficient
              typename RepresenterType::DatasetPointerType  sampleShape = m_asmodel->GetStatisticalModel()->DrawSample(currentSample.GetCoefficients());
              typename RepresenterType::DatasetPointerType sample = m_asmodel->GetRepresenter()->TransformMesh(sampleShape, currentSample.GetRigidTransform());


            for( int i = 0; i < m_targetPoints.size(); ++i) {

              PointType closestPtOnSample =  m_closestPoint->findClosestPoint(sample, m_targetPoints[i]);
              double d = m_eval->evalSample(closestPtOnSample-m_targetPoints[i]);

                distance += d;
            }
            return distance;

          }
        private:
          const vector< PointType > m_targetPoints;
          const ActiveShapeModelType* m_asmodel;
          PositionEvaluator* m_eval;
          const ClosestPointType* m_closestPoint;

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




      template <class T>
      class ASMEvaluator : public DistributionEvaluator< ChainSampleType > {
        typedef  ASMPreprocessedImage<typename RepresenterType::DatasetType> PreprocessedImageType;

      public:
          ASMEvaluator(ActiveShapeModelType* asmodel, PreprocessedImageType* image, PointSamplerType* sampler) :
                  m_asmodel(asmodel),
                  m_image(image),
                  m_sampler(sampler)
          {}

          // DistributionEvaluatorInterface interface
      public:
          virtual double evalSample(const ChainSampleType& currentSample) {



              typename RepresenterType::DatasetPointerType  sampleShape = m_asmodel->GetStatisticalModel()->DrawSample(currentSample.GetCoefficients());
              typename RepresenterType::DatasetPointerType currentModelInstance = m_asmodel->GetRepresenter()->TransformMesh(sampleShape, currentSample.GetRigidTransform());


              unsigned numProfilePointsUsed = 300;
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


          ASMLikelihoodForChunk evalSampleForProfiles(unsigned profileIdStart, unsigned profileIdEnd, unsigned step, typename RepresenterType::DatasetPointerType currentModelInstance, const ChainSampleType& currentSample) {

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


      // TODO: Full image evaluator, i.e. ASM Evaluator is missing



      template <class T>
      class HURangeASMEvaluator : public DistributionEvaluator< ChainSampleType > {

      public:
        HURangeASMEvaluator(ActiveShapeModelType* asmodel, PreprocessedImageType* image, const FeatureExtractorType* imageValueExtractor) :
          m_asmodel(asmodel),
          m_image(image),
          m_extractor(imageValueExtractor) {}

        // DistributionEvaluatorInterface interface
      public:
        virtual double evalSample(const ChainSampleType& currentSample) {

          
          typename RepresenterType::DatasetPointerType  sampleShape = m_asmodel->GetStatisticalModel()->DrawSample(currentSample.GetCoefficients());
          typename RepresenterType::DatasetPointerType currentModelInstance = m_asmodel->GetRepresenter()->TransformMesh(sampleShape, currentSample.GetRigidTransform());

          
          std::cout << "Evaluating current sample in HU Evaluator." << std::endl;


          unsigned numProfilePointsUsed = 300;
          unsigned step = m_asmodel->GetProfiles().size() / numProfilePointsUsed;
          unsigned numChunks = boost::thread::hardware_concurrency() + 1;
          std::vector<boost::future<ASMLikelihoodForChunk>* > futvec;


          for (unsigned i = 0; i < numChunks; ++i) {
            unsigned profileIdStart = m_asmodel->GetProfiles().size() / numChunks  * i;
            unsigned profileIdEnd = m_asmodel->GetProfiles().size() / numChunks * (i + 1);

            boost::future<ASMLikelihoodForChunk> *fut = new boost::future<ASMLikelihoodForChunk>(
              boost::async(boost::launch::async, boost::bind(&HURangeASMEvaluator<T>::evalSampleForProfiles,
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


        ASMLikelihoodForChunk evalSampleForProfiles(unsigned profileIdStart, unsigned profileIdEnd, unsigned step, 
          typename RepresenterType::DatasetPointerType currentModelInstance, const ChainSampleType& currentSample) {

          
          FeatureExtractorType* fe = m_extractor->CloneForTarget(m_asmodel, currentSample.GetCoefficients(), currentSample.GetRigidTransform());

          double loglikelihood = 0.0;

          for (unsigned i = profileIdStart; i < profileIdEnd; i += step) {

            ASMProfile profile = m_asmodel->GetProfiles()[i];
            long ptId = profile.GetPointId();

            statismo::VectorType features;


            bool ok = fe->ExtractFeatures(features, m_image, currentModelInstance->GetPoint(ptId));

            if (ok) {
              const int intercept = -1024;
              const int threshold = 400;
              // Tests have shown at least for the Esophagus example, that the 
              // first half of the points are inside.
              for (size_t j = 0; j < features.size() / 2; j++) {
                if (features[j] + intercept > threshold)
                {
                  std::cout << "Value above threshold: " << features[j] + intercept << std::endl;
                  auto diff = features[j] + intercept - threshold;
                  loglikelihood += diff*diff;
                }
              }

            }
            else {
              std::cout << "feature not ok " << std::endl;
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
        const FeatureExtractorType* m_extractor;
        RandomGenerator* m_rGen;
      };



    public:
      // "Script"
      /*============================================================================*/
      template <class T>
      class BasicSampling {
        public:
          static MarkovChain<ChainSampleType >* buildChain(
                  const Representer<T>* representer,
                  const ClosestPoint<typename Representer<T>::DatasetPointerType, typename Representer<T>::PointType>* closestPoint,
                  MHFittingConfiguration config,
                  ASMPreprocessedImage<typename Representer<T>::DatasetType>* preprocessedFeatureImage,
                  ASMPreprocessedImage<typename Representer<T>::DatasetType>* interpolatedImage,
                  const FeatureExtractorType* imageValueExtractor,
              vector<PointType> targetPoints,
              ActiveShapeModelType* asmodel,
                  PointSamplerType* asmPointSampler,
              RigidTransformPointerType transform,
              statismo::VectorType coeffs) {

            // basics
            RandomGenerator* rGen = new RandomGenerator(42);
            ChainSampleType init = ChainSampleType(coeffs,transform);

            // Proposals // TODO: add ASMProposal
            GaussianModelUpdate* poseProposalRough = new GaussianModelUpdate(0.2,rGen);
            GaussianModelUpdate* poseProposalFine = new GaussianModelUpdate(0.05,rGen);
              GaussianModelUpdate* poseProposalFinest = new GaussianModelUpdate(0.01,rGen);
              ASMModelUpdate* asmProposal = new ASMModelUpdate(config.GetAsmFittingconfiguration(), asmodel, preprocessedFeatureImage, asmPointSampler);


            vector< typename RandomProposal< ChainSampleType >::GeneratorPair> gaussMixtureProposalVector(3);
            gaussMixtureProposalVector[0] = pair<ProposalGenerator<ChainSampleType >*,double>(poseProposalRough,0.2);
            gaussMixtureProposalVector[1] = pair<ProposalGenerator<ChainSampleType >*,double>(poseProposalFine,0.2);
              gaussMixtureProposalVector[2] = pair<ProposalGenerator<ChainSampleType >*,double>(poseProposalFinest,0.6);

              RandomProposal<ChainSampleType >* gaussMixtureProposal = new RandomProposal<ChainSampleType >(gaussMixtureProposalVector,rGen);

            Gaussian3DPositionDifferenceEvaluator* diffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 1.0);
            PointEvaluator<T>* pointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints,asmodel,diffEval);
            ModelPriorEvaluator* modelPriorEvaluator = new ModelPriorEvaluator();
//
            auto HUEvaluator = new HURangeASMEvaluator<T>(asmodel, interpolatedImage, imageValueExtractor);

              QuietLogger <ChainSampleType>* ql = new QuietLogger<ChainSampleType>();
              MarkovChain<ChainSampleType >* landMarkchain =  new MetropolisHastings<ChainSampleType >(gaussMixtureProposal, pointEval, ql, init, rGen );
              MarkovChainProposal<ChainSampleType>* lmChainProposal = new MarkovChainProposal<ChainSampleType>(landMarkchain, 10);


//            Gaussian3DPositionDifferenceEvaluator* wideDiffEval = new Gaussian3DPositionDifferenceEvaluator(asmodel->GetRepresenter(), 5.0);
//            PointEvaluator<T>* widePointEval = new PointEvaluator<T>(representer, closestPoint, targetPoints,asmodel,wideDiffEval);
            ASMEvaluator<T>* asmEvaluator = new ASMEvaluator<T>(asmodel, preprocessedFeatureImage, asmPointSampler);

              std::vector<DistributionEvaluator<ChainSampleType >*> evaluatorList;
              vector< typename RandomProposal< ChainSampleType >::GeneratorPair> finalProposalVector;


              if (targetPoints.size() > 0) {
                  evaluatorList.push_back(pointEval);
                  evaluatorList.push_back(asmEvaluator);
                  evaluatorList.push_back(modelPriorEvaluator);
                  evaluatorList.push_back(HUEvaluator);

                  finalProposalVector.push_back(pair<ProposalGenerator<ChainSampleType >*,double>(lmChainProposal,0.9));
                  finalProposalVector.push_back(pair<ProposalGenerator<ChainSampleType >*,double>(asmProposal,0.1));

              } else {
                  evaluatorList.push_back(asmEvaluator);
                  evaluatorList.push_back(modelPriorEvaluator);
                  evaluatorList.push_back(HUEvaluator);

                  finalProposalVector.push_back(pair<ProposalGenerator<ChainSampleType >*,double>(gaussMixtureProposal,1.0));
                  finalProposalVector.push_back(pair<ProposalGenerator<ChainSampleType >*,double>(asmProposal,0.0));

              }
              RandomProposal<ChainSampleType >* finalProposal = new RandomProposal<ChainSampleType >(finalProposalVector,rGen);



              ProductEvaluator<ChainSampleType >* productEvaluator = new ProductEvaluator<ChainSampleType >(evaluatorList);

            // markov chain
              BestMatchLogger<ChainSampleType>* bml = new BestMatchLogger<ChainSampleType>();

              MarkovChain<ChainSampleType >* chain = new MetropolisHastings<ChainSampleType>(finalProposal, productEvaluator, ql, init, rGen);

            return chain;
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
      typedef MHFittingResult<RigidTransformPointerType> ChainSampleType;

      MHFittingStep(MarkovChain<ChainSampleType >* chain)
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


      static MHFittingStep *Create(MarkovChain<ChainSampleType >* chain) {
        return new MHFittingStep( chain);
      }

      ~MHFittingStep() {
      }

      ChainSampleType Perform() const {


        //          ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
        //          ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();


        // runs a markov chain and returns only the accepted proposals
        ChainSampleType params;//(m_sourceTransform,m_sourceCoefficients);
        m_chain->next(params);
        ChainSampleType mhResult(params.GetCoefficients(), params.GetRigidTransform());
        return mhResult;
      }

    private:
      MarkovChain<ChainSampleType >* m_chain;
  };
}

#endif //STATISMO_ASMFITTING_H

