//
// Created by luetma00 on 15.11.16.
//

#ifndef STATISMO_EVALUATORS_H_H
#define STATISMO_EVALUATORS_H_H

#include <stdexcept>

#include "MHFittingParameters.h"
#include "MeshOperations.h"
#include "Sampling/DistributionEvaluator.h"
#include "Types.h"

namespace mhfitting {


    using sampling::DistributionEvaluator;


    // Evaluators
    /*============================================================================*/

    class PositionEvaluator : public DistributionEvaluator<PointType>
    {
    public:
        /// dummy gradient, zero, do nothing
        virtual void evalGradient( PointType& gradient, PointType const& sample )
        {
            throw std::logic_error("no gradient implemented in PositionEvaluator");
        }

        /// Returns whether this evaluator is separable into x and y coordinates
        virtual bool isSeparable() const
        {
            return false;
        }
    };


    template <class PointType>
    class ModelPriorEvaluator : public DistributionEvaluator<MHFittingParameters>
    {
    public:
        virtual double evalSample(const MHFittingParameters& currentSample) {

            statismo::VectorType coefficients = currentSample.GetCoefficients();

            return -0.5 * coefficients.size() *log( 2*M_PI ) - 0.5 * coefficients.squaredNorm();

        }

        /// dummy gradient, zero, do nothing
        virtual void evalGradient( PointType& gradient, PointType const& sample )
        {
            throw std::logic_error("no gradient implemented in PositionEvaluator");
        }

        /// Returns whether this evaluator is separable into x and y coordinates
        virtual bool isSeparable() const
        {
            return false;
        }
    };


    template <class T>
    class InLungOrBoneEvaluator : public DistributionEvaluator< MHFittingParameters > {
        typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> ClosestPointType;


        typedef statismo::Representer<T> RepresenterType;
        typedef statismo::StatisticalModel<T> StatisticalModelType;
        typedef typename RepresenterType::PointType PointType;

    public:
        InLungOrBoneEvaluator( const statismo::Representer<T>* representer, const ClosestPointType* closestPoint, const StatisticalModelType* model) :
                m_model(model),
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
            unsigned numPoints = m_model->GetRepresenter()->GetDomain().GetNumberOfPoints();
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
        const std::vector< PointType > m_targetPoints;
        const StatisticalModelType* m_model;
        PositionEvaluator* m_eval;
        const ClosestPointType* m_meshOperations;

    };



    template <class T>
    class Gaussian3DPositionDifferenceEvaluator : public PositionEvaluator
    {
    private:
        typedef statismo::Representer<T> RepresenterType;
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
        typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> ClosestPointType;

    public:
        PointEvaluator( const statismo::Representer<T>* representer, const ClosestPointType* closestPoint, const CorrespondencePoints& correspondencePoints, ActiveShapeModelType* asmodel, PositionEvaluator* evaluator) :
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
        typedef MeshOperations<typename statismo::Representer<T>::DatasetPointerType, typename statismo::Representer<T>::PointType> ClosestPointType;
        typedef std::map<float, std::vector<statismo::VectorType> > LineMapType;
    public:
        LineEvaluator( const statismo::Representer<T>* representer, const ClosestPointType* closestPoint, const std::vector< PointType >& targetPoints, ActiveShapeModelType* asmodel, double noisevar) :
                m_targetPoints(targetPoints),
                m_asmodel(asmodel),
                m_closestPoint(closestPoint),
                m_noiseVar(noisevar)
        {

            // TODO HACK, we separate the points into lines
            for (unsigned i = 0 ; i < targetPoints.size(); ++i ) {
                statismo::VectorType pt = representer->PointToVector(targetPoints[i]);
                if (m_lineMap.find(pt(2)) != m_lineMap.end()) {
                    m_lineMap.insert(std::make_pair(pt(2), std::vector<statismo::VectorType>()));
                }
                m_lineMap[pt(2)].push_back(pt);

            }
            std::cout << "number of lines " << m_lineMap.size();
            statismo::VectorType mean = statismo::VectorType::Zero(1);
            mean(0) = 0.0;
            statismo::MatrixType cov = statismo::MatrixType::Zero(1,1);
            cov(0,0)=m_noiseVar;
            m_likelihoodModel = statismo::MultiVariateNormalDistribution(mean, cov);
        }

        // DistributionEvaluatorInterface interface
    public:


        double likelihoodForLine(const std::vector<statismo::VectorType> & line, typename RepresenterType::DatasetPointerType sample) {
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
                statismo::VectorType distAsVec(1);
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
        statismo::MultiVariateNormalDistribution m_likelihoodModel;
        const std::vector< PointType > m_targetPoints;
        const ActiveShapeModelType* m_asmodel;
        const ClosestPointType* m_closestPoint;
        LineMapType m_lineMap;
        double m_noiseVar;
    };


}

#endif //STATISMO_EVALUATORS_H_H
