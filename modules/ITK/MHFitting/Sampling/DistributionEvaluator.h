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
 *         @file DistributionEvaluator.h
 *
 *        @brief Standard implementation of the DistributionEvaluator
 *
 *         @date 5 Sep 2012
 *      @authors Sandro Schönborn (ses)\n
 *               sandro.schoenborn@unibas.ch\n
 *               University of Basel, Switzerland
 */
/*============================================================================*/

#ifndef DISTRIBUTIONEVALUATOR_H
#define DISTRIBUTIONEVALUATOR_H

#include <string>
#include <stdexcept>

#include "MarkovChainInterfaces.h"

namespace sampling
{
  /** \brief Distribution Evaluator - for Markov Chain - Interface description -- ATTENTION: always LOG probabilities! */
  template <typename T>
  class DistributionEvaluator :
      public virtual DistributionEvaluatorInterface<T>,
      public virtual GradientEvaluatorInterface<T>
  {
    public:
      typedef T SampleType;

      DistributionEvaluator() :
        name_("AnonymousDistributionEvaluator")
      {}

      DistributionEvaluator( std::string const& name ) :
        name_( name )
      {}

      virtual ~DistributionEvaluator() {}

      /// The standard implementation does not have a gradient implemented
      virtual void evalGradient(SampleType &gradient, const SampleType &sample)
      {
        throw std::runtime_error("No gradient implemented in DistributionEvaluator");
      }

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

#endif // DISTRIBUTIONEVALUATOR_H
