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
/*
 * =====================================================================================
 *
 *       Filename:  ForwardChain.h
 *
 *    Description:  Implementation of a simple Foward MarkovChain
 *
 *        Version:  1.0
 *        Created:  13.04.2010 15:56:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */

#ifndef  FORWARDCHAIN_INC
#define  FORWARDCHAIN_INC

#include "../MarkovChain.h"

namespace sampling
{

  /** \brief Implementation of a simple forward Markov Chain
   *
   * Proposals are generated by the ProposalGenerator and are evaluated
   * by the DistributionEvaluator.
   * Template T: type of samples (can be abstract)
   * */
  template <typename T>
  class ForwardChain : public MarkovChain<T>
  {
    public:
      typedef T SampleType;

    public:
      typedef ProposalGenerator<T> Generator;
      typedef DistributionEvaluator<T> Evaluator;
      typedef ChainLogger<T> Logger;

    public:
      ForwardChain( Generator* propgen, Evaluator* disteval, Logger* chainlogger, const T& initsample ) :
        currentSample( initsample ),
        pGenerator( propgen ),
        pEvaluator( disteval ),
        pLogger( chainlogger )
      {
        currentProbVal = pEvaluator->evalSample( currentSample );
        pLogger->notifyReset( currentSample, currentProbVal, pGenerator, pEvaluator );
      }

      virtual ~ForwardChain() {}

      void current( SampleType& rSample ) const
      {
        rSample = currentSample;
      }

      double currentValue() const
      {
        return currentProbVal;
      }

      /** \brief Get next sample: Sample from Proposal Distribution and Evaluate
      */
      void next( SampleType& rNextSample )
      {
        // generate a Proposal
        pGenerator->generateProposal( rNextSample, currentSample );

        // Evaluate Proposal
        double probVal = pEvaluator->evalSample( rNextSample );

        currentSample = rNextSample;
        currentProbVal = probVal;
        pLogger->notifyAccept( currentSample, probVal, pGenerator, pEvaluator  );
      }

      double setState( const SampleType& newState )
      {
        currentSample = newState;
        currentProbVal = pEvaluator->evalSample( currentSample );
        pLogger->notifyReset( currentSample, currentProbVal, pGenerator, pEvaluator  );
        return currentProbVal;
      }

    private:
      SampleType currentSample;
      double currentProbVal;

      Generator* pGenerator;
      Evaluator* pEvaluator;
      Logger* pLogger;
  };
}

#endif   /* ----- #ifndef FORWARDCHAIN_INC  ----- */
