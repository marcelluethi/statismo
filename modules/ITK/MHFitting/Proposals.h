//
// Created by luetma00 on 16.11.16.
//

#ifndef STATISMO_PROPOSALS_H
#define STATISMO_PROPOSALS_H

#include <stdexcept>
#include <memory>

#include "MHFittingParameters.h"
#include "MeshOperations.h"
#include "Sampling/ProposalGenerator.h"
#include "Sampling/Proposal/RandomProposal.h"
#include "Sampling/ProposalGenerator.h"
#include "Sampling/RandomGenerator.h"
#include "Types.h"

namespace mhfitting {


    using sampling::ProposalGenerator;
    using sampling::RandomGenerator;
    using statismo::VectorType;
    using sampling::RandomProposal;


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


            proposal = MHFittingParameters(newShapeParams, currentSample.GetRigidTransformParameters(), currentSample.GetRotationCenter());
        }

        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }

    private:
        double sigmaShape;
        RandomGenerator* rgen;
        unsigned m_maxNumberOfShapeParameters;
        double sigmaTranslation;
        double sigmaRotation;

    };


    class RandomShapeProposal : public ProposalGenerator<MHFittingParameters > {
        typedef RandomProposal< MHFittingParameters > RandomProposalType;

    public:
        RandomShapeProposal( unsigned numPCAComponents, RandomGenerator* rgen )
                : m_numPCAComponents(numPCAComponents), m_rGen(rgen) {

            RandomShapeUpdate* shapeUpdateRough = new RandomShapeUpdate(0.1, numPCAComponents, m_rGen);
            RandomShapeUpdate* shapeUpdateFine = new RandomShapeUpdate(0.05, numPCAComponents, m_rGen);
            RandomShapeUpdate* shapeUpdateFinest = new RandomShapeUpdate(0.01, numPCAComponents, m_rGen);


            std::vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> gaussMixtureProposalVector;
            gaussMixtureProposalVector.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateRough, 0.1));
            gaussMixtureProposalVector.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateFine, 0.3));
            gaussMixtureProposalVector.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(shapeUpdateFinest, 0.6));
            m_mixtureProposal.reset(new RandomProposalType(gaussMixtureProposalVector, m_rGen));
        }

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){
            VectorType shapeParams = currentSample.GetCoefficients();

            m_mixtureProposal->generateProposal(proposal, currentSample);
        }

        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }

    private:
        unsigned m_numPCAComponents;
        RandomGenerator* m_rGen;
        std::unique_ptr<RandomProposalType> m_mixtureProposal;

    };



    class RotationUpdate : public ProposalGenerator<MHFittingParameters > {
    private:
        unsigned m_axis;
        double m_sigmaRotation;
        RandomGenerator* m_rgen;

    public:
        RotationUpdate( unsigned axis,double stepSizeRotation, RandomGenerator* rgen)
                : m_axis(axis) , m_sigmaRotation(stepSizeRotation),  m_rgen(rgen) {}

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            VectorType newRigidParams = currentSample.GetRigidTransformParameters();
            newRigidParams[m_axis] += m_rgen->normalDbl() * m_sigmaRotation;
            proposal = MHFittingParameters( currentSample.GetCoefficients(), newRigidParams, currentSample.GetRotationCenter());
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
        TranslationUpdate(double stepSizeTranslations, unsigned axis, RandomGenerator* rgen)
                : m_axis(axis) , m_sigmaTranslation(stepSizeTranslations),  m_rgen(rgen) {}

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            VectorType newRigidParams = currentSample.GetRigidTransformParameters();
            newRigidParams[3 + m_axis] += m_rgen->normalDbl() * m_sigmaTranslation;
            proposal = MHFittingParameters( currentSample.GetCoefficients(), newRigidParams, currentSample.GetRotationCenter());
        }

        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }
    };

    class RandomPoseProposal : public ProposalGenerator<MHFittingParameters> {

    public:
        RandomPoseProposal( RandomGenerator* rgen) : m_rgen(rgen) {



            RotationUpdate* rotUpdateX = new RotationUpdate(0, 0.01, m_rgen);
            RotationUpdate* rotUpdateY = new RotationUpdate(1, 0.01, m_rgen);
            RotationUpdate* rotUpdateZ = new RotationUpdate(2, 0.01, m_rgen);

            std::vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> rotUpdateVec;
            rotUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateX, 0.2));
            rotUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateY, 0.2));
            rotUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdateZ, 0.6));
            RandomProposal<MHFittingParameters >* rotUpdate = new RandomProposal<MHFittingParameters >(rotUpdateVec, m_rgen);


            TranslationUpdate* transUpdateXSmall = new TranslationUpdate(0, 1, m_rgen);
            TranslationUpdate* transUpdateYSmall = new TranslationUpdate(1, 1, m_rgen);
            TranslationUpdate* transUpdateZSmall = new TranslationUpdate(2, 1, m_rgen);
            TranslationUpdate* transUpdateXLarge = new TranslationUpdate(0, 10, m_rgen);
            TranslationUpdate* transUpdateYLarge = new TranslationUpdate(1, 10, m_rgen);
            TranslationUpdate* transUpdateZLarge = new TranslationUpdate(2, 10, m_rgen);


            std::vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> transUpdateVec;
            transUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateXSmall, 0.1));
            transUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateYSmall, 0.1));
            transUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateZSmall, 0.1));
            transUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateXLarge, 0.1));
            transUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateYLarge, 0.1));
            transUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdateZLarge, 0.1));

            RandomProposal<MHFittingParameters >* transUpdate = new RandomProposal<MHFittingParameters >(transUpdateVec, m_rgen);

            std::vector< typename RandomProposal< MHFittingParameters >::GeneratorPair> poseUpdateVec;
            poseUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdate, 0.1));
            poseUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdate, 0.1));

            m_randomPoseProposal.reset(new RandomProposal<MHFittingParameters >(poseUpdateVec, m_rgen));

        }

        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample) {

//                    val pitchProposalC = RotationUpdate(0.05f, 1)
//                    val pitchProposalI = RotationUpdate(0.01f, 1)
//                    val pitchProposalF = RotationUpdate(0.001f, 1)
//                    val rotationPitch = MixtureProposal(0.1 *: pitchProposalC + 0.4 *: pitchProposalI + 0.5 *: pitchProposalF)
//
//                    val rollProposalC = RotationUpdate(0.2f, 0)
//                    val rollProposalI = RotationUpdate(0.1f, 0)
//                    val rollProposalF = RotationUpdate(0.01f, 0)
//                    val rotationRoll = MixtureProposal(0.1 *: rollProposalC + 0.4 *: rollProposalI + 0.5 *: rollProposalF)
//
//
//                    val rotationProposal = MixtureProposal(0.5 *: rotationRoll + 0.2 *: rotationPitch + 0.2 *: rotationYaw)
//
//
//                    val translationXProposalC = TranslationUpdate(10f, 0)
//                    val translationXProposalI = TranslationUpdate(5f, 0)
//                    val translationXProposalF = TranslationUpdate(1f, 0)
//                    val translationXProposal = MixtureProposal(0.1 *: translationXProposalC + 0.4 *: translationXProposalI + 0.5 *: translationXProposalF)
//
//
//                    val translationYProposalC = TranslationUpdate(10f, 1)
//                    val translationYProposalI = TranslationUpdate(5f, 1)
//                    val translationYProposalF = TranslationUpdate(1f, 1)
//                    val translationYProposal = MixtureProposal(0.1 *: translationYProposalC + 0.4 *: translationYProposalI + 0.5 *: translationYProposalF)
//
//
//                    val translationZProposalC = TranslationUpdate(10f, 2)
//                    val translationZProposalI = TranslationUpdate(5f, 2)
//                    val translationZProposalF = TranslationUpdate(1f, 2)
//                    val translationZProposal = MixtureProposal(0.1 *: translationZProposalC + 0.4 *: translationZProposalI + 0.5 *: translationZProposalF)
//
//
//                    val translationProposal = MixtureProposal(0.3 *: translationXProposal + 0.3 *: translationYProposal + 0.4 *: translationZProposal)
//
//
//                    val poseProposal = MixtureProposal(0.4 *: rotationProposal + 0.4 *: translationProposal)
//
//                    poseProposal
            m_randomPoseProposal->generateProposal(proposal, currentSample);

        }


        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end) {
            return 0.0;
        }

    private:
        RandomGenerator* m_rgen;
        std::unique_ptr<RandomProposal<MHFittingParameters > > m_randomPoseProposal;
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






    template <class T>
    class RigidICPProposal : public ProposalGenerator<MHFittingParameters > {
        typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> MeshOperationsType;
        typedef std::vector<std::pair<PointType, PointType> > CorrespondencePointsType;
        typedef std::vector<PointType> LinePointsType;

    public:
        RigidICPProposal(const MeshOperationsType* meshOperations, const CorrespondencePointsType correspondencePoints, const LinePointsType linePoints) : m_meshOperations(meshOperations) {


            for (CorrespondencePointsType::const_iterator it = correspondencePoints.begin(); it != m_correspondencePoints.end(); ++it) {
                m_targetPoints.push_back(it->second);
            }
            for (LinePointsType::const_iterator it = linePoints.begin(); it != linePoints.end(); ++it) {
                m_targetPoints.push_back(*it);
            }
        }

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            proposal = m_meshOperations->rigidICP(currentSample, m_targetPoints);
        }


        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }

    private:
        const MeshOperationsType* m_meshOperations;
        CorrespondencePointsType m_correspondencePoints;
        LinePointsType m_targetPoints;
    };



} // namespace mhfitting

#endif //STATISMO_PROPOSALS_H
