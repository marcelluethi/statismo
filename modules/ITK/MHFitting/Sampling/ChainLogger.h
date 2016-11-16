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
 *         @file ChainLogger.h
 *
 *        @brief Standard implementation of a ChainLogger
 *
 *         @date 5 Sep 2012
 *      @authors Sandro Schönborn (ses)\n
 *               sandro.schoenborn@unibas.ch\n
 *               University of Basel, Switzerland
 */
/*============================================================================*/

#ifndef CHAINLOGGER_H
#define CHAINLOGGER_H

#include <string>
#include "MarkovChainInterfaces.h"

namespace sampling
{
  /** \brief Chain logger facility - Interface description with empty/quiet standard implementation */
  template <typename T>
  class ChainLogger : public virtual ChainLoggerInterface<T>
  {
    public:
      typedef T SampleType;
      typedef ProposalGeneratorInterface<T>     Generator;
      typedef DistributionEvaluatorInterface<T> Evaluator;
      typedef ChainLoggerInterface<T>           Logger;

      ChainLogger() :
        name_("AnonymousChainLogger")
      {}

      ChainLogger( std::string const& name ) :
        name_(name)
      {}

      virtual ~ChainLogger() {}

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
#endif // CHAINLOGGER_H
