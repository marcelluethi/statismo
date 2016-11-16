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
 *         @file ProposalGenerator.h
 *
 *        @brief Standard implementation of a ProposalGenerator
 *
 *         @date 5 Sep 2012
 *      @authors Sandro Schönborn (ses)\n
 *               sandro.schoenborn@unibas.ch\n
 *               University of Basel, Switzerland
 */
/*============================================================================*/

#ifndef PROPOSALGENERATOR_H
#define PROPOSALGENERATOR_H

#include <string>
#include "MarkovChainInterfaces.h"

namespace sampling
{
  /** \brief Proposal Generator for Markov Chain processes - Standard implementation */
  template <typename T>
  class ProposalGenerator :
      public virtual ProposalGeneratorInterface<T>
  {
    public:
      typedef T SampleType;

    public:
      ProposalGenerator() :
        name_("AnonymousProposalGenerator")
      {}

      ProposalGenerator( std::string const& name ) :
        name_(name)
      {}

      virtual ~ProposalGenerator() {}

      virtual std::string getName() const
      {
        return name_;
      }

      virtual void setName( std::string const& name )
      {
        name_ = name;
      }

      virtual double transitionRatio(const SampleType &start, const SampleType &end)
      {
        return this->transitionProbability(start,end) - this->transitionProbability(end,start);
      }

    private:
      std::string name_;
  };
}
#endif // PROPOSALGENERATOR_H
