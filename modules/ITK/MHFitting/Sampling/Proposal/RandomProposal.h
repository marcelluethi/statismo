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
 *       Filename:  RandomProposal.h
 *
 *    Description:  Choose a Proposal at random, multinomial distribution
 *
 *        Version:  1.0
 *        Created:  02/25/2011 09:01:52 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */


#ifndef  RANDOMPROPOSAL_INC
#define  RANDOMPROPOSAL_INC

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

#include "../RandomGenerator.h"
#include "../MarkovChain.h"

namespace sampling
{
  /** \brief MixtureProposal: Generates a proposal by choosing a random sub generator,
    * the sub generator is guaranteed to change only on "generateProposal" or "add"
    */
  template <typename T>
  class RandomProposal : public ProposalGenerator<T>
  {
    public:
      typedef ProposalGenerator<T> Generator;
      typedef typename std::vector<Generator*>::iterator Iterator;
      typedef std::pair<Generator*,double> GeneratorPair;

    public:
      RandomProposal( std::vector<GeneratorPair> generators, RandomGenerator* rndgenerator ) :
        vecGens(generators),
        pRndEngine( rndgenerator )
      {
        normalize();
        itCurrent = 0;
      }

      void generateProposal( T& rSample, T const& stCurrent )
      {
        advance();
        vecGens[itCurrent].first->generateProposal( rSample, stCurrent );
      }

      double transitionProbability( T const& from, T const& to )
      {
        // Calculate transitions for whole mixture (does not need current generator)
        // log sum for mixture distribution
        std::vector<double> logProbs( vecGens.size() );
        for( size_t i = 0; i < vecGens.size(); ++i )
          logProbs[i] = vecGens[i].first->transitionProbability(from,to);

        // for the log sum biggest exponent use: log(mixturevalue*exp(logProbs)) = log(mixturevalue)+logProbs
        std::vector<double> logExps( logProbs );
        for (size_t i = 0; i < vecGens.size(); ++i)
          logExps[i] = logProbs[i] + std::log(vecGens[i].second);

        double maxExp = std::max( *(std::max_element(logExps.begin(), logExps.end())), -1e300 ); // no -inf for maxExp

        double totalProbability = 0.0;
        for ( size_t i = 0; i < vecGens.size(); ++i )
          totalProbability += vecGens[i].second * std::exp( logProbs[i] - maxExp );

        return log(totalProbability) + maxExp;
      }

      /// Returns the name of the actual proposal generator
      virtual std::string getName() const
      {
        if ( vecGens.size() > 0 )
          return vecGens.at(itCurrent).first->getName();
        else
          return ProposalGenerator<T>::getName();
      }

    public:
      std::vector<GeneratorPair> const& getGenerators() const
      {
        return vecGens;
      }

    private:
      /// Normalize the probabilities to sum up to 1
      void normalize()
      {
        double dSum = 0.0;
        for( size_t i = 0; i < vecGens.size(); i++ )
          dSum += vecGens[i].second;
        for( size_t i = 0; i < vecGens.size(); i++ )
          vecGens[i].second = vecGens[i].second / dSum;
      }

      /// Advance to the next Proposal Generator - automatically called by generateProposal
      void advance()
      {
        if ( vecGens.size() == 0 )
        {
          throw std::runtime_error("Random Proposals: No Generators registered!");
        }

        double dRnd = pRndEngine->uniformDbl();
        double dCSum = 0.0;
        itCurrent = 0;

        for( itCurrent = 0; itCurrent < vecGens.size() && dRnd > dCSum; itCurrent++ )
          dCSum += vecGens[itCurrent].second;
        itCurrent--;
      }

      std::vector<GeneratorPair> vecGens;
      RandomGenerator* pRndEngine;
      size_t itCurrent;
  };
}
#endif   /* ----- #ifndef RANDOMPROPOSAL_INC  ----- */
