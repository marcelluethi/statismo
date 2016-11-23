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

#include "Sampling/DistributionEvaluator.h"
#include "Sampling/Evaluator/ProductEvaluator.h"
#include "Sampling/MarkovChain.h"
#include "Sampling/Algorithm/ForwardChain.h"
#include "Sampling/Algorithm/Metropolis.h"
#include "Sampling/Logger/QuietLogger.h"
#include "Sampling/Algorithm/MetropolisHastings.h"
#include "Sampling/Logger/BestMatchLogger.h"
#include "Sampling/Proposal/MarkovChainProposal.h"


#include "Evaluators.h"
#include "Proposals.h"

namespace mhfitting {

    using namespace sampling;
    using namespace std;

    using statismo::VectorType;
    using statismo::MatrixType;
    using statismo::MultiVariateNormalDistribution;
    using statismo::Representer;


    class MHFittingConfiguration {
    private:

    public:
        MHFittingConfiguration(const statismo::ASMFittingConfiguration &asmFittingConfiguration,
                               unsigned numberOfProfilePoints = 100) : m_asmFittingConfiguration(
                asmFittingConfiguration) {}

        const statismo::ASMFittingConfiguration &
        GetAsmFittingconfiguration() const { return m_asmFittingConfiguration; }

    private:
        statismo::ASMFittingConfiguration m_asmFittingConfiguration;
    };



//
//    template<class T>
//    class UncertaintyLogger : public ChainLogger<MHFittingParameters> {
//        typedef unsigned int PointId;
//
//    public:
//        typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> MeshOperationsType;
//
//            UncertaintyLogger(const statismo::Representer<T> *representer, const MeshOperationsType *meshOps,
//                              const CorrespondencePoints &correspondencePoints, const vector<PointType> &targetPoints)
//
//                    : m_representer(representer),
//                      m_meshOps(meshOps),
//                      m_correspondencePoints(correspondencePoints),
//                      m_targetPoints(targetPoints) {
//                // fill the map with an empty vector of points
//                for (unsigned i = 0; i < m_correspondencePoints.size(); ++i) {
//                    unsigned id = m_correspondencePoints[i].first;
//                    m_uncertaintyPerCorrespondingPoint.insert(std::make_pair(id, std::vector<PointType>()));
//                }
//
//
//            }
//
//            ~UncertaintyLogger() {
//                std::cout << "destructor of uncertainty logger" << std::endl;
//
//            }
//
//
//            virtual void notifyAccept(
//                    const MHFittingParameters &parameters,
//                    const double &dProbValue,
//                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
//                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
//            ) {
//
//                typename statismo::Representer<T>::DatasetPointerType sample = m_meshOps->transformMesh(parameters);
//                for (unsigned i = 0; i < m_correspondencePoints.size(); ++i) {
//                    unsigned id = m_correspondencePoints[i].first;
//                    PointType pointOnSample = m_meshOps->getPointWithId(sample, id);
//                    m_uncertaintyPerCorrespondingPoint[id].push_back(pointOnSample);
//                }
//
//                // we keep the parameters, such that we can compute the average at the end
//                m_params.push_back(parameters);
//            }
//
//            /** \brief Function gets called whenever an algorithm rejects a proposal */
//            virtual void notifyReject(
//                    const MHFittingParameters &sample,
//                    const double &dProbValue,
//                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
//                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
//            ) {}
//
//            /** \brief Function gets called whenever an algorithm is reset to a new state */
//            virtual void notifyReset(
//                    const MHFittingParameters &state,
//                    const double &dProbValue,
//                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
//                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
//            ) {}
//
//
//            statismo::MultiVariateNormalDistribution uncertaintyForCorrespondencePoint(unsigned id) {
//
//                if (m_params.size() == 0) {
//                    throw std::runtime_error("we did not have an accepted sample in the unceratinty estimation");
//                }
//                statismo::VectorType avgRotationParameters = statismo::VectorType::Zero(
//                        m_params.front().GetRotationParameters().GetParameters().rows());
//                statismo::VectorType avgTranslationParameters = statismo::VectorType::Zero(
//                        m_params.front().GetTranslationParameters().GetParameters().rows());
//                for (typename std::list<MHFittingParameters>::const_iterator it = m_params.begin();
//                     it != m_params.end(); ++it) {
//                    avgRotationParameters += it->GetRotationParameters().GetParameters();
//                    avgTranslationParameters += it->GetTranslationParameters().GetParameters();
//                }
//                avgRotationParameters /= m_params.size();
//                avgTranslationParameters /= m_params.size();
//
//                statismo::VectorType mean = statismo::VectorType::Zero(3);
//                statismo::MatrixType cov = statismo::MatrixType::Zero(3, 3);
//
//
//                const std::vector<PointType> &v = m_uncertaintyPerCorrespondingPoint[id];
//
//
//                for (typename std::vector<PointType>::const_iterator it = v.begin(); it != v.end(); ++it) {
//                    statismo::VectorType ptAsVec = m_representer->PointToVector(
//                            m_meshOps->transformToModelSpace(avgRotationParameters, avgTranslationParameters, *it));
//                    mean += ptAsVec;
//                }
//                mean /= v.size();
//
//
//                for (typename std::vector<PointType>::const_iterator it = v.begin(); it != v.end(); ++it) {
//                    VectorType ptAsVec = m_representer->PointToVector(
//                            m_meshOps->transformToModelSpace(avgRotationParameters, avgTranslationParameters, *it));
//                    cov += (ptAsVec - mean) * (ptAsVec - mean).transpose();
//                }
//                cov /= (v.size() - 1);
//
//                MultiVariateNormalDistribution mvn(mean, cov);
//                return mvn;
//
//            }
//
//        private:
//            const statismo::Representer<T> *m_representer;
//            const MeshOperationsType *m_meshOps;
//            std::vector<PointType> m_targetPoints;
//            CorrespondencePoints m_correspondencePoints;
//            std::list<MHFittingParameters> m_params;
//            std::map<PointId, std::vector<PointType> > m_uncertaintyPerCorrespondingPoint;
//        };


        template<class T>
        class InLungLogger : public ChainLogger<MHFittingParameters> {
            typedef typename statismo::Representer<T>::PointType PointType;
            typedef unsigned int PointId;

        public:


            typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> MeshOperationsType;

            InLungLogger(ChainLogger<MHFittingParameters> *otherLogger, const statismo::Representer<T> *representer,
                         const MeshOperationsType *meshOps, const char *loggerName)

                    : m_compositeLogger(otherLogger),
                      m_representer(representer),
                      m_meshOps(meshOps),
                      m_loggerName(loggerName) {
            }

            ~InLungLogger() {
                std::cout << "destructor of uncertainty logger" << std::endl;

            }


            virtual void notifyAccept(
                    const MHFittingParameters &parameters,
                    const double &dProbValue,
                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
            ) {


                m_compositeLogger->notifyAccept(parameters, dProbValue, proposal, evaluator);


                typename Representer<T>::DatasetPointerType sample = m_meshOps->transformMesh(parameters);
                unsigned numPoints = m_representer->GetDomain().GetNumberOfPoints();
                bool isInsideBone = false;
                unsigned numInBone = 0;
                for (int i = 0; i < numPoints; ++i) {

                    if (m_meshOps->huAtPoint(sample, i) > 1150 || m_meshOps->huAtPoint(sample, i) < 400) {
//                      isInsideBone = true;
//                      break;
                        numInBone += 1;
                    }
                }

//              if (isInsideBone) {
//                  distance = -std::numeric_limits<double>::infinity();
//
//              } else {
//                  distance = 0;
//              }

                std::cout << "Accepted Step in " << m_loggerName << " : New number in bone is " << numInBone
                          << " parameter norm " << parameters.GetCoefficients().GetParameters().squaredNorm() << std::endl;

            }

            /** \brief Function gets called whenever an algorithm rejects a proposal */
            virtual void notifyReject(
                    const MHFittingParameters &sample,
                    const double &dProbValue,
                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
            ) {
                m_compositeLogger->notifyReject(sample, dProbValue, proposal, evaluator);

            }

            /** \brief Function gets called whenever an algorithm is reset to a new state */
            virtual void notifyReset(
                    const MHFittingParameters &state,
                    const double &dProbValue,
                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
            ) {
                m_compositeLogger->notifyReset(state, dProbValue, proposal, evaluator);
            }


        private:
            ChainLogger<MHFittingParameters> *m_compositeLogger;
            const statismo::Representer<T> *m_representer;
            const MeshOperationsType *m_meshOps;
            const char *m_loggerName;
        };


        template<class T>
        class PoseLogger : public ChainLogger<MHFittingParameters> {

        public:
            PoseLogger(const char *loggerName)
                    : m_loggerName(loggerName) {
            }


            virtual void notifyAccept(
                    const MHFittingParameters &parameters,
                    const double &dProbValue,
                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
            ) {


                std::cout << "accepted propseal in logger " << m_loggerName << std::endl;


            }

            /** \brief Function gets called whenever an algorithm rejects a proposal */
            virtual void notifyReject(
                    const MHFittingParameters &sample,
                    const double &dProbValue,
                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
            ) {

                std::cout << "rejected propseal in logger " << m_loggerName << std::endl;
            }

            /** \brief Function gets called whenever an algorithm is reset to a new state */
            virtual void notifyReset(
                    const MHFittingParameters &state,
                    const double &dProbValue,
                    ProposalGeneratorInterface<MHFittingParameters> *proposal,
                    DistributionEvaluatorInterface<MHFittingParameters> *evaluator
            ) {}


        private:
            const char *m_loggerName;
        };


        struct ASMLikelihoodForChunk {


            ASMLikelihoodForChunk(double _aggregatedLikelihood) : aggregatedLikelihood(_aggregatedLikelihood) {}

            double aggregatedLikelihood;

            // emulate move semantics, as boost::async seems to depend on it.
            ASMLikelihoodForChunk &operator=(BOOST_COPY_ASSIGN_REF(ASMLikelihoodForChunk)rhs) { // Copy assignment

                if (&rhs != this) {
                    copyMembers(rhs);
                }
                return *this;
            }

            ASMLikelihoodForChunk(BOOST_RV_REF(ASMLikelihoodForChunk)that) { //Move constructor
                copyMembers(that);
            }

            ASMLikelihoodForChunk &operator=(BOOST_RV_REF(ASMLikelihoodForChunk)rhs) { //Move assignment
                if (&rhs != this) {
                    copyMembers(rhs);
                }
                return *this;
            }

        private:
            BOOST_COPYABLE_AND_MOVABLE(ASMLikelihoodForChunk)

            void copyMembers(const ASMLikelihoodForChunk &that) {
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



        // "Script"
        /*============================================================================*/
        template<class T>
        class BasicSampling {
        public:


//            // estimate the uncertainty at the given correspondence Points, by sampling from the initialPoseChain.
//            static std::map<unsigned, statismo::MultiVariateNormalDistribution>
//            estimatePointUncertaintyForInitialPoseChain(const statismo::Representer<T> *representer,
//                                                        const MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> *meshOperations,
//                                                        CorrespondencePoints correspondencePoints,
//                                                        vector<PointType> targetPoints,
//                                                        ActiveShapeModelType *asmodel,
//                                                        MHFittingParameters &initialParameters) {
///                UncertaintyLogger<T> *ul = new UncertaintyLogger<T>(representer, meshOperations, correspondencePoints,
//                                                                    targetPoints);
//                MarkovChain <MHFittingParameters> *chain = buildInitialPoseChain(representer, meshOperations,
//                                                                                 correspondencePoints, targetPoints,
//                                                                                 asmodel, initialParameters, ul);
//
//                MHFittingParameters params;//(m_sourceTransform,m_sourceCoefficients);
//                for (unsigned i = 0; i < 1000; ++i) {
//                    chain->next(params);
//                }
//
//                std::map<unsigned, MultiVariateNormalDistribution> uncertaintyMap;
//
//                for (unsigned i = 0; i < correspondencePoints.size(); ++i) {
//                    unsigned id = correspondencePoints[i].first;
//
//                    uncertaintyMap.insert(std::make_pair(id, ul->uncertaintyForCorrespondencePoint(id)));
//                }
//                return uncertaintyMap;
//
//            }

            static MarkovChain<MHFittingParameters> *buildInitialPoseChain(
                    const statismo::Representer<T> *representer,
                    const MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> *meshOperations,
                    CorrespondencePoints correspondencePoints,
                    vector<PointType> targetPoints,
                    ActiveShapeModelType *asmodel,
                    MHFittingParameters &initialParameters,
                    ChainLogger<MHFittingParameters> *chainLogger) {


                // basics
                RandomGenerator *rGen = new RandomGenerator(42);
                RandomPoseProposal *randomPoseProposal = new RandomPoseProposal(rGen);
                randomPoseProposal->setName("RandomPoseProposal");
                RigidICPProposal<MeshType>* rigidICPProposal = new RigidICPProposal<MeshType>(meshOperations, correspondencePoints, targetPoints);
                rigidICPProposal->setName("rigidICPProposal");

                std::vector<typename RandomProposal<MHFittingParameters>::GeneratorPair> poseMixtureProposalVec;
                poseMixtureProposalVec.push_back(std::make_pair(randomPoseProposal, 0.5));
                poseMixtureProposalVec.push_back(std::make_pair(rigidICPProposal, 0.5));
                RandomProposal<MHFittingParameters>* poseMixtureProposal = new RandomProposal<MHFittingParameters>(poseMixtureProposalVec, rGen);

                Gaussian3DPositionDifferenceEvaluator<MeshType> *diffEval = new Gaussian3DPositionDifferenceEvaluator<MeshType>(asmodel->GetRepresenter(), 2.0);
                PointEvaluator<T> *pointEval = new PointEvaluator<T>(representer, meshOperations, correspondencePoints,
                                                                     asmodel, diffEval);
                LineEvaluator<T> *lineEval = new LineEvaluator<T>(representer, meshOperations, targetPoints, asmodel, 1.0);

                ModelPriorEvaluator<MeshType> *modelPriorEvaluator = new ModelPriorEvaluator<MeshType>();


                std::vector<DistributionEvaluator<MHFittingParameters> *> evaluatorList;
                evaluatorList.push_back(pointEval);
                evaluatorList.push_back(lineEval);
                evaluatorList.push_back(modelPriorEvaluator);

                MarkovChain <MHFittingParameters> *lmChain = new MetropolisHastings<MHFittingParameters>(
                        poseMixtureProposal, new ProductEvaluator<MHFittingParameters>(evaluatorList), chainLogger,
                        initialParameters, rGen);

                return lmChain;
            }



            static MarkovChain<MHFittingParameters> *buildPoseAndShapeChain(
                    const statismo::Representer<T> *representer,
                    const MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> *closestPoint,
                    CorrespondencePoints correspondencePoints,
                    vector<PointType> targetPoints,
                    ActiveShapeModelType *asmodel,
                    MHFittingParameters &initialParameters,
                    ChainLogger<MHFittingParameters> *logger) {


                // basics
                RandomGenerator *rGen = new RandomGenerator(42);
                RandomPoseProposal *randomPoseProposal = new RandomPoseProposal(rGen);
                RandomShapeProposal *randomShapeProposal = new RandomShapeProposal(
                        asmodel->GetStatisticalModel()->GetNumberOfPrincipalComponents(), rGen);
                std::vector<typename RandomProposal<MHFittingParameters>::GeneratorPair> poseAndShapeMixtureVec;
                poseAndShapeMixtureVec.push_back(std::make_pair(randomPoseProposal, 0.5));
                poseAndShapeMixtureVec.push_back(std::make_pair(randomShapeProposal, 0.5));
                RandomProposal<MHFittingParameters> *poseAndShapeProposal = new RandomProposal<MHFittingParameters>(
                        poseAndShapeMixtureVec, rGen);

                Gaussian3DPositionDifferenceEvaluator<MeshType> *diffEval = new Gaussian3DPositionDifferenceEvaluator<MeshType>(
                        asmodel->GetRepresenter(), 1.0);
                PointEvaluator<T> *pointEval = new PointEvaluator<T>(representer, closestPoint, correspondencePoints,
                                                                     asmodel, diffEval);
                LineEvaluator<T> *lineEval = new LineEvaluator<T>(representer, closestPoint, targetPoints, asmodel,
                                                                  3.0);

                ModelPriorEvaluator<MeshType> *modelPriorEvaluator = new ModelPriorEvaluator<MeshType>();

                std::vector<DistributionEvaluator<MHFittingParameters> *> evaluatorList;
                evaluatorList.push_back(pointEval);
                evaluatorList.push_back(lineEval);
                evaluatorList.push_back(modelPriorEvaluator);

                MarkovChain <MHFittingParameters> *lmChain = new MetropolisHastings<MHFittingParameters>(
                        poseAndShapeProposal, new ProductEvaluator<MHFittingParameters>(evaluatorList), logger,
                        initialParameters, rGen);

                return lmChain;
            }


            static MarkovChain<MHFittingParameters> *buildLmAndHuChain(
                    const statismo::Representer<T> *representer,
                    const MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> *meshOperations,
                    CorrespondencePoints correspondencePoints,
                    vector<PointType> targetPoints,
                    ActiveShapeModelType *asmodel,
                    MHFittingParameters &initialParameters,
                    ChainLogger<MHFittingParameters> *logger) {

                unsigned numPCAComponents = asmodel->GetStatisticalModel()->GetNumberOfPrincipalComponents();

                // basics
                RandomGenerator *rGen = new RandomGenerator(42);
                MHFittingParameters init = MHFittingParameters(initialParameters.GetCoefficients(),
                                                               initialParameters.GetRotationParameters(),
                                                               initialParameters.GetTranslationParameters(),
                                                               initialParameters.GetRotationCenter());

                RandomShapeUpdate *shapeUpdateRough = new RandomShapeUpdate(0.2, numPCAComponents, rGen);
                RandomShapeUpdate *shapeUpdateFine = new RandomShapeUpdate(0.1, numPCAComponents, rGen);
                RandomShapeUpdate *shapeUpdateFinest = new RandomShapeUpdate(0.025, numPCAComponents, rGen);
                RotationUpdate *rotUpdateX = new RotationUpdate(0, 0.01, rGen);
                RotationUpdate *rotUpdateY = new RotationUpdate(1, 0.01, rGen);
                RotationUpdate *rotUpdateZ = new RotationUpdate(2, 0.01, rGen);

                vector<typename RandomProposal<MHFittingParameters>::GeneratorPair> rotUpdateVec(3);
                rotUpdateVec[0] = pair<ProposalGenerator<MHFittingParameters> *, double>(rotUpdateX, 0.8);
                rotUpdateVec[1] = pair<ProposalGenerator<MHFittingParameters> *, double>(rotUpdateY, 0.1);
                rotUpdateVec[2] = pair<ProposalGenerator<MHFittingParameters> *, double>(rotUpdateZ, 0.1);
                RandomProposal<MHFittingParameters> *rotUpdate = new RandomProposal<MHFittingParameters>(rotUpdateVec,
                                                                                                         rGen);


                vector<typename RandomProposal<MHFittingParameters>::GeneratorPair> gaussMixtureProposalVector(4);
                gaussMixtureProposalVector[0] = pair<ProposalGenerator<MHFittingParameters> *, double>(shapeUpdateRough,
                                                                                                       0.1);
                gaussMixtureProposalVector[1] = pair<ProposalGenerator<MHFittingParameters> *, double>(shapeUpdateFine,
                                                                                                       0.2);
                gaussMixtureProposalVector[2] = pair<ProposalGenerator<MHFittingParameters> *, double>(
                        shapeUpdateFinest,
                        0.4);
                gaussMixtureProposalVector[3] = pair<ProposalGenerator<MHFittingParameters> *, double>(rotUpdate, 0.2);


                RandomProposal<MHFittingParameters> *gaussMixtureProposal = new RandomProposal<MHFittingParameters>(
                        gaussMixtureProposalVector, rGen);

                Gaussian3DPositionDifferenceEvaluator<MeshType> *diffEval = new Gaussian3DPositionDifferenceEvaluator<MeshType>(
                        asmodel->GetRepresenter(), 1.0);
                PointEvaluator<T> *pointEval = new PointEvaluator<T>(representer, meshOperations, correspondencePoints,
                                                                     asmodel, diffEval);
                LineEvaluator<T> *lineEval = new LineEvaluator<T>(representer, meshOperations, targetPoints, asmodel,
                                                                  0.1);

                ModelPriorEvaluator<PointType> *modelPriorEvaluator = new ModelPriorEvaluator<PointType>();

                InLungOrBoneEvaluator<T> *huEvaluator = new InLungOrBoneEvaluator<T>(representer, meshOperations,
                                                                                     asmodel->GetStatisticalModel());

                std::vector<DistributionEvaluator<MHFittingParameters> *> huEvaluatorList;
                //  huEvaluatorList.push_back(huEvaluator);
                //huEvaluatorList.push_back(modelPriorEvaluator);
//              huEvaluatorList.push_back(pointEval);
                huEvaluatorList.push_back(lineEval);


                QuietLogger <MHFittingParameters> *ql = new QuietLogger<MHFittingParameters>();
                InLungLogger<T> *loggerFilterChain = new InLungLogger<T>(logger, representer, meshOperations,
                                                                         "filter chain");
                MarkovChain <MHFittingParameters> *huChain = new MetropolisHastings<MHFittingParameters>(
                        gaussMixtureProposal, new ProductEvaluator<MHFittingParameters>(huEvaluatorList), ql, init,
                        rGen);
                MarkovChainProposal <MHFittingParameters> *huChainProposal = new MarkovChainProposal<MHFittingParameters>(
                        huChain, 1, true);


                std::vector<DistributionEvaluator<MHFittingParameters> *> lmAndHuEvaluatorList;
                //lmAndHuEvaluatorList.push_back(pointEval);
                //lmAndHuEvaluatorList.push_back(lineEval);
                lmAndHuEvaluatorList.push_back(huEvaluator);
                lmAndHuEvaluatorList.push_back(modelPriorEvaluator);

                //InLungLogger <T>* loggerFinalChain = new InLungLogger<T>(representer, meshOperations, "final chain");
                MarkovChain <MHFittingParameters> *lmAndHuChain = new MetropolisHastings<MHFittingParameters>(
                        huChainProposal, new ProductEvaluator<MHFittingParameters>(lmAndHuEvaluatorList),
                        loggerFilterChain, init, rGen);
                return lmAndHuChain;
            }


        };


        template<typename TPointSet, typename TImage>
        class MHFittingStep {
            typedef statismo::ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
            typedef statismo::ASMPointSampler<TPointSet, TImage> PointSamplerType;
            typedef typename statismo::Representer<TPointSet>::PointType PointType;
            typedef statismo::StatisticalModel<TPointSet> StatisticalModelType;
            typedef statismo::ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
            typedef statismo::ASMPreprocessedImage<TPointSet> PreprocessedImageType;
            typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
            typedef statismo::ASMFittingStep<TPointSet, TImage> ASMFittingStepType;

            MHFittingStep(MarkovChain<MHFittingParameters> *chain)
                    :
                    m_chain(chain) {}

        private:
            class ProfileResult {
            public:
                unsigned int pointId;
                PointType candidatePoint;
                PointType transformedCandidatePoint;
            };


        public:


            static MHFittingStep *Create(MarkovChain<MHFittingParameters> *chain) {
                return new MHFittingStep(chain);
            }

            ~MHFittingStep() {
            }

            MHFittingParameters Perform() const {


                //          ASMFittingStepType* asmFittingStep = ASMFittingStepType::Create(m_configuration.GetAsmFittingconfiguration(), m_model, m_sourceCoefficients, m_sourceTransform, m_target, m_sampler);
                //          ASMFittingResult<RigidTransformPointerType> result = asmFittingStep->Perform();


                // runs a markov chain and returns only the accepted proposals
                MHFittingParameters params;//(m_sourceTransform,m_sourceCoefficients);
                m_chain->next(params);
                MHFittingParameters mhResult(params.GetCoefficients(), params.GetRotationParameters(), params.GetTranslationParameters(),
                                             params.GetRotationCenter());
                return mhResult;
            }

        private:
            MarkovChain<MHFittingParameters> *m_chain;


        };


}

#endif //STATISMO_ASMFITTING_H


