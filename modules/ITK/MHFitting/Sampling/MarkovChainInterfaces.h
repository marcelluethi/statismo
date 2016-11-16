/*
 * Copyright 2015 University of Basel, Graphics and Vision Research Group
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*============================================================================*/
/**
 *         @file MarkovChainInterfaces.h
 *
 *        @brief \todo
 *
 *         @date 5 Sep 2012
 *      @authors Sandro Schönborn (ses)\n
 *               sandro.schoenborn@unibas.ch\n
 *               University of Basel, Switzerland
 */
/*============================================================================*/

#ifndef MARKOVCHAININTERFACES_H
#define MARKOVCHAININTERFACES_H

#include "NamedObjectInterface.h"

// \TODO make methods const where possible!

namespace sampling
{
  /** \brief Proposal Generator for Markov Chain processes - Interface */
  template <typename T>
  class ProposalGeneratorInterface : public virtual NamedObjectInterface
  {
    public:
      typedef T SampleType;

      virtual ~ProposalGeneratorInterface() {}

      /** \brief Function generating a proposal given the current chain state */
      virtual void generateProposal(
        SampleType& rProposal,            //< [out] New proposed state value
        const SampleType& currentSample   //< [in] Current State (Propsal might depend on this)
      ) = 0;

      /** \brief Function returning the transition probability between states */
      virtual double transitionProbability(
        SampleType const& start,  //< [in] state to start from
        SampleType const& end     //< [out] state to end at
      ) = 0;

      /** \brief Function returning the transition ratio between states FBR */
      virtual double transitionRatio(
        SampleType const& start,  //< [in] state to start from
        SampleType const& end     //< [out] state to end at
      ) = 0;
  };

  /** \brief Distribution Evaluator - for Markov Chain - Interface description -- ATTENTION: always LOG probabilities! */
  template <typename T>
  class DistributionEvaluatorInterface : public virtual NamedObjectInterface
  {
    public:
      typedef T SampleType;

      virtual ~DistributionEvaluatorInterface() {}

      /** \brief Return the probability value of this sample */
      virtual double evalSample( const SampleType& currentSample ) = 0;

  };

  /** \brief Gradient Evaluator calculates gradients */
  template <typename T>
  class GradientEvaluatorInterface : public virtual NamedObjectInterface
  {
    public:
      typedef T SampleType;

      virtual ~GradientEvaluatorInterface() {}

      /** \brief Calculate the gradient at the current sample */
      virtual void evalGradient(
          SampleType& gradient,       /**< [OUT] holds the gradient after the call */
          SampleType const& sample    /**< [IN] the sample to evaluate the gradient at */
      ) = 0;
  };

  /** \brief Chain logger facility - Interface description with empty/quiet standard implementation */
  template <typename T>
  class ChainLoggerInterface : public virtual NamedObjectInterface
  {
    public:
      typedef T SampleType;

      virtual ~ChainLoggerInterface() {}

      /** \brief Function gets called whenever an algorithm accepts a proposal */
      virtual void notifyAccept(
          const T& sample,
          const double& dProbValue,
          ProposalGeneratorInterface<T>* proposal,
          DistributionEvaluatorInterface<T>* evaluator
      ) = 0;

      /** \brief Function gets called whenever an algorithm rejects a proposal */
      virtual void notifyReject(
          const T& sample,
          const double& dProbValue,
          ProposalGeneratorInterface<T>* proposal,
          DistributionEvaluatorInterface<T>* evaluator
          ) = 0;

      /** \brief Function gets called whenever an algorithm is reset to a new state */
      virtual void notifyReset (
          const T& state,
          const double& dProbValue,
          ProposalGeneratorInterface<T>* proposal,
          DistributionEvaluatorInterface<T>* evaluator
          ) = 0;
  };

  /** \brief A Scalable interface, needed for some Proposals */
  class ScalableInterface : public virtual NamedObjectInterface
  {
    public:
      virtual void scale( double factor ) = 0;
  };

  /** \brief Very simple Markov Chain interface */
  template <typename T>
  class MarkovChainInterface : public virtual NamedObjectInterface
  {
    public:
      /** \brief The data type of the samples */
      typedef T SampleType;

      typedef ProposalGeneratorInterface<T> Generator;
      typedef DistributionEvaluatorInterface<T> Evaluator;
      typedef ChainLoggerInterface<T> Logger;

    public:
      virtual ~MarkovChainInterface() {}

      /** \brief Get next sample from chain */
      virtual void next(
        SampleType& rNextSample /**< [out] Holds next sample after call */
      ) = 0;

      /** \brief Get current sample, do not advance */
      virtual void current(
          SampleType& rSample /**< [out] Holds the current sample after call */
      ) const = 0;

      /** \brief Return Probability value of current sample - depends on implementation, can be a log-prob */
      virtual double currentValue() const = 0;

      /** \brief Set the current state of the Chain */
      virtual double setState(
          const SampleType& newState /**< [IN] The state the chain is set to */
      ) = 0;
  };

}
#endif // MARKOVCHAININTERFACES_H
