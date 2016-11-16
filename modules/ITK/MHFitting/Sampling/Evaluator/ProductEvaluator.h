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
 *       Filename:  ProductEvaluator.h
 *
 *    Description:  Evaluate a Product
 *
 *        Version:  1.0
 *        Created:  02/01/2011 08:57:19 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */

#ifndef  PRODUCTEVALUATOR_INC
#define  PRODUCTEVALUATOR_INC

#include <vector>
#include <algorithm>
#include <numeric>

#include "../MarkovChain.h"

namespace sampling
{
  /** \brief Log Version of Product of terms, e.g. Likelihood * Prior */
  template <typename T>
  class ProductEvaluator : public DistributionEvaluator<T>
  {
    public:
      typedef DistributionEvaluator<T> Evaluator;
      typedef std::vector<Evaluator*> EvaluatorList;
      typedef typename EvaluatorList::iterator EvaluatorIterator;

    public:
      ProductEvaluator(EvaluatorList evaluators) :
        Evaluators(evaluators)
      {}

      virtual ~ProductEvaluator() {}

      double evalSample( const T& current )
      {
        double sum = 0.0;
        for( size_t i = 0; i < Evaluators.size(); i++ )
        {
          double dS = Evaluators[i]->evalSample( current );
          sum += dS;
        }
        return sum;
      }

      void evalGradient( T& gradient, T const& current )
      {
        for( EvaluatorIterator itE = Evaluators.begin(); itE != Evaluators.end(); itE++ )
        {
          T gradCurrent( gradient );
          (*itE)->evalGradient( gradCurrent, current );

          for( size_t i = 0; i < gradient.size(); ++i )
            gradient[i] += gradCurrent[i];
        }
      }

      virtual std::string getName() const
      {
        if ( Evaluators.size() == 1 )
          return Evaluators.at(0)->getName();

        if ( Evaluators.size() > 0 )
        {
          //std::string fullName = DistributionEvaluator<T>::getName() + ":" + Evaluators.at(0)->getName();
          std::string fullName = Evaluators.at(0)->getName();
          for ( size_t i = 1; i < Evaluators.size(); ++i )
            fullName += std::string("+")+Evaluators[i]->getName();
          return fullName;
        }
        return DistributionEvaluator<T>::getName();
      }

    protected:
      std::vector<Evaluator*> Evaluators;
  };

  /** \brief Log Version of Product of terms, e.g. Likelihood * Prior */
  template <typename T>
  class TrimmedProductEvaluator : public ProductEvaluator<T>
  {
    public:
      typedef typename ProductEvaluator<T>::EvaluatorIterator EvaluatorIterator;
      typedef typename ProductEvaluator<T>::EvaluatorList EvaluatorList;

    public:
      TrimmedProductEvaluator(EvaluatorList evaluators, double useFraction = 1.0) :
        ProductEvaluator<T>(evaluators),
        dUseFraction(useFraction)
      {}

      virtual ~TrimmedProductEvaluator() {}

      double evalSample( const T& current )
      {
        double sum = 0.0;
        int iNEvaluators = this->Evaluators.size();
        int dropNEvaluators = (int)round((1.0-dUseFraction)*iNEvaluators);

        // fill vector with values
        std::vector<double> dS(iNEvaluators);
        for( EvaluatorIterator itE = this->Evaluators.begin(); itE != this->Evaluators.end(); itE++ )
        {
          dS.push_back((*itE)->evalSample( current ));
        }
        // sort vector and sum only a part
        std::sort(dS.begin(),dS.end());
        sum += std::accumulate(dS.begin()+dropNEvaluators,dS.end(),0.0);

//        double dropsum  = std::accumulate(dS.begin(), dS.begin()+dropNEvaluators,0.0);
//        double totalsum = std::accumulate(dS.begin(), dS.end(), 0.0);
//        std::cout << "TrimmedProductEvaluator: used sum=" << sum << ", ignored=" << dropsum << ", total=" << totalsum << "\n";
        return sum;
      }

      // TODO evalGradient

    private:
      double dUseFraction;
  };

  /** \brief Log Version of Product of terms, e.g. Likelihood * Prior */
  template <typename T>
  class LimitedProductEvaluator : public ProductEvaluator<T>
  {
    public:
      typedef typename ProductEvaluator<T>::EvaluatorIterator EvaluatorIterator;
      typedef typename ProductEvaluator<T>::EvaluatorList EvaluatorList;

    public:
      LimitedProductEvaluator(EvaluatorList evaluators, double limit) :
        ProductEvaluator<T>(evaluators),
        dLimit(limit)
      {}

      virtual ~LimitedProductEvaluator() {}

      double evalSample( const T& current )
      {
        double sum = 0.0;

        // fill vector with values
        for( EvaluatorIterator itE = this->Evaluators.begin(); itE != this->Evaluators.end(); itE++ )
        {
          double dVal = (*itE)->evalSample( current );
          sum += std::max<double>( dLimit, dVal );
          //std::cout << "LimitedProductEvaluator: eval value=" << dVal << ", limited: " << std::max<double>( dLimit, dVal ) << "\n";
        }

        return sum;
      }

      // TODO evalGradient

    private:
      double dLimit;
  };
}

#endif   /* ----- #ifndef PRODUCTEVALUATOR_INC  ----- */

