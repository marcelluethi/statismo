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
 *       Filename:  BestMatchLogger.h
 *
 *    Description:  Stores the highest probability sample so far
 *
 *        Version:  1.0
 *        Created:  02/01/2011 09:33:23 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */


#ifndef  BESTMATCHLOGGER_INC
#define  BESTMATCHLOGGER_INC

#include "../MarkovChain.h"

#include <cmath>
#include <limits>

namespace sampling
{
  /** \brief Keep best sample so far */
  template <typename T>
  class BestMatchLogger : public ChainLogger<T>
  {
    public:
      typedef typename ChainLogger<T>::SampleType SampleType;
      typedef typename ChainLogger<T>::Generator  Generator;
      typedef typename ChainLogger<T>::Evaluator  Evaluator;
      typedef typename ChainLogger<T>::Logger     Logger;

    public:
      BestMatchLogger() :
        best(),
        dBestProbVal(0)
      {}

      // Interface methods
      virtual void notifyAccept( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
      {
        if ( std::isfinite(pValue) && pValue > dBestProbVal )
        {
          best = cparm;
          dBestProbVal = pValue;
        }
      }

      virtual void notifyReject( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
      {
        // Empty so far
      }

      /// Reset state
      virtual void notifyReset( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
      {
        if ( std::isfinite(pValue) ) {
          best = cparm;
          dBestProbVal = pValue;
        }
        else { // keep invalid value as -inifity, will be replaced by next valid sample
          best = cparm;
          dBestProbVal = -std::numeric_limits<double>::infinity();
        }
      }

      /// Returns a std::copy of the best sample
      T getBest() const
      {
        return best;
      }

      /// Returns the probability value of the best sample
      double getBestProbability() const
      {
        return dBestProbVal;
      }

    private:
      T best;
      double dBestProbVal;
  };
}
#endif   /* ----- #ifndef BESTMATCHLOGGER_INC  ----- */
