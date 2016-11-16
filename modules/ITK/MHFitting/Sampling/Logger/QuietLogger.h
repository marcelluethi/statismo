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
 *       Filename:  QuietLogger.h
 *
 *    Description:  Quiet Logger, does not produce any output
 *
 *        Version:  1.0
 *        Created:  02/01/2011 09:31:33 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */


#ifndef  QUIETLOGGER_INC
#define  QUIETLOGGER_INC

#include "../MarkovChain.h"

namespace sampling
{
  /** \brief Empty Logger, does nothing */
  template <typename T>
  class QuietLogger : public ChainLogger<T>
  {
    public:
      typedef typename ChainLogger<T>::SampleType SampleType;
      typedef typename ChainLogger<T>::Generator  Generator;
      typedef typename ChainLogger<T>::Evaluator  Evaluator;
      typedef typename ChainLogger<T>::Logger     Logger;

    public:
      QuietLogger() {}
      virtual ~QuietLogger() {}

      void notifyAccept( const T& rSample, const double& pValue, Generator* proposal, Evaluator* evaluator ) {}
      void notifyReject( const T& rSample, const double& pValue, Generator* proposal, Evaluator* evaluator ) {}
      void notifyReset ( const T& rSample, const double& pValue, Generator* proposal, Evaluator* evaluator ) {}
  };
}
#endif   /* ----- #ifndef QUIETLOGGER_INC  ----- */
