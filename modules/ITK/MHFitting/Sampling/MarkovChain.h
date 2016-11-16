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
 *       Filename:  MarkovChain.h
 *
 *    Description:  A very simple Markov Chain framework for sampling applications
 *
 *        Version:  1.0
 *        Created:  13.04.2010 15:59:48
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */

#ifndef MARKOVCHAIN_H
#define MARKOVCHAIN_H

#include "NamedObjectInterface.h"
#include "MarkovChainInterfaces.h"

#include "ProposalGenerator.h"
#include "DistributionEvaluator.h"
#include "ChainLogger.h"

namespace sampling
{
  /** \brief Very simple Markov Chain interface */
  template <typename T>
  class MarkovChain : public virtual MarkovChainInterface<T>
  {
    public:
      MarkovChain() :
        name_("AnonymousMarkovChain")
      {}

      MarkovChain( std::string const& name ) :
        name_(name)
      {}

      virtual ~MarkovChain() {}

      virtual std::string getName() const
      {
        return name_;
      }

      virtual void setName( std::string const& name )
      {
        name_ = name;
      }

    private:
      std::string name_;
  };
}

#endif
