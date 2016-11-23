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
            VectorType shapeParams = currentSample.GetCoefficients().GetParameters();
            VectorType newShapeParams(shapeParams.size());
            newShapeParams.fill(0);
            for (unsigned i = 0; i < std::min(m_maxNumberOfShapeParameters, static_cast<unsigned>(shapeParams.size())); ++i) {
                newShapeParams[i] = shapeParams[i] + rgen->normalDbl() * sigmaShape;
            }


            proposal = MHFittingParameters(newShapeParams, currentSample.GetRotationParameters(), currentSample.GetTranslationParameters(), currentSample.GetRotationCenter());
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

            this->setName("randomShapeProposal");

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
            VectorType shapeParams = currentSample.GetCoefficients().GetParameters();

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

            VectorType newRotationParams = currentSample.GetRotationParameters().GetParameters();
            newRotationParams[m_axis] += m_rgen->normalDbl() * m_sigmaRotation;
            proposal = MHFittingParameters( currentSample.GetCoefficients(), newRotationParams, currentSample.GetTranslationParameters(), currentSample.GetRotationCenter());
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
        TranslationUpdate( unsigned axis, double stepSizeTranslations, RandomGenerator* rgen)
                : m_axis(axis) , m_sigmaTranslation(stepSizeTranslations),  m_rgen(rgen) {}

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            VectorType newTranslationParameters = currentSample.GetTranslationParameters().GetParameters();
            newTranslationParameters[m_axis] += m_rgen->normalDbl() * m_sigmaTranslation;
            proposal = MHFittingParameters( currentSample.GetCoefficients(), currentSample.GetRotationParameters(), newTranslationParameters, currentSample.GetRotationCenter());
        }

        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }
    };

    class RandomPoseProposal : public ProposalGenerator<MHFittingParameters> {

    public:
        RandomPoseProposal( RandomGenerator* rgen) : m_rgen(rgen) {

            this->setName("randomPoseProposal");


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
            poseUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(rotUpdate, 0.5));
            poseUpdateVec.push_back(std::pair<ProposalGenerator<MHFittingParameters >*, double>(transUpdate, 0.5));

            m_randomPoseProposal.reset(new RandomProposal<MHFittingParameters >(poseUpdateVec, m_rgen));


        }

        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample) {

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
      public:                rigidICPProposal->setName("rigidICPProposal");
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

        typedef std::vector<PointType> LinePoints;

    public:
        RigidICPProposal(const MeshOperationsType* meshOperations, const CorrespondencePoints& correspondencePoints, const LinePoints& linePoints, double stepLength) : m_meshOperations(meshOperations), m_stepLength(stepLength) {

            this->setName("RigidICPProposal");
            for (CorrespondencePoints::const_iterator it = correspondencePoints.begin(); it != correspondencePoints.end(); ++it) {
                m_targetPoints.push_back(it->second);
            }
            for (LinePoints::const_iterator it = linePoints.begin(); it != linePoints.end(); ++it) {
                m_targetPoints.push_back(*it);
            }
        }

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            MHFittingParameters::RigidParameters newProposal = m_meshOperations->alignModelInstanceToTargetPoints(currentSample, m_targetPoints);

            statismo::VectorType newRotationParameter = currentSample.GetRotationParameters().GetParameters() + (newProposal.GetRotationParameters().GetParameters() - currentSample.GetRotationParameters().GetParameters()) * m_stepLength;
            statismo::VectorType newTranslationParameter = currentSample.GetTranslationParameters().GetParameters() + (newProposal.GetTranslationParameters().GetParameters() - currentSample.GetTranslationParameters().GetParameters()) * m_stepLength;
            proposal = MHFittingParameters(currentSample.GetCoefficients(),
                                           MHFittingParameters::RotationParameters(newRotationParameter),
                                           MHFittingParameters::TranslationParameters(newTranslationParameter),
                                           currentSample.GetRotationCenter());

        }


        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }

    private:
        const MeshOperationsType* m_meshOperations;
        LinePoints m_targetPoints;
        double m_stepLength;
    };




    template <class T>
    class ShapeICPProposal : public ProposalGenerator<MHFittingParameters > {
        typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> MeshOperationsType;

        typedef std::vector<PointType> LinePoints;

    public:
        ShapeICPProposal(const MeshOperationsType* meshOperations, const CorrespondencePoints& correspondencePoints, const LinePoints& linePoints, double stepLength) :
                m_meshOperations(meshOperations),
                m_correspondenPoints(correspondencePoints),
                m_linePoints(linePoints),
                m_stepLength(stepLength) {

            this->setName("ShapeICPProposal");
        }

        // ProposalGeneratorInterface interface
    public:
        virtual void generateProposal(MHFittingParameters& proposal, const MHFittingParameters& currentSample){

            MHFittingParameters::ModelParameters newProposal = m_meshOperations->projectModelInstanceToTargetPoints(currentSample, m_correspondenPoints, m_linePoints);

            statismo::VectorType newModelParameters = currentSample.GetCoefficients().GetParameters() + (newProposal.GetParameters() - currentSample.GetCoefficients().GetParameters()) * m_stepLength;

            proposal = MHFittingParameters(newModelParameters,
                                           currentSample.GetRotationParameters(),
                                           currentSample.GetTranslationParameters(),
                                           currentSample.GetRotationCenter());

        }


        virtual double transitionProbability(const MHFittingParameters& start, const MHFittingParameters& end){
            return 0.0;
        }

    private:
        const MeshOperationsType* m_meshOperations;
        CorrespondencePoints m_correspondenPoints;
        LinePoints m_linePoints;
        double m_stepLength;
    };




} // namespace mhfitting

#endif //STATISMO_PROPOSALS_H
