/*
 * Copyright 2016 University of Basel, Switzerland
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
 *
 */

#ifndef  LOGGERDISTRIBUTOR_INC
#define  LOGGERDISTRIBUTOR_INC

#include <vector>

#include "../MarkovChain.h"

namespace sampling
{
    /** \brief Copy logging stream to multiple loggers - Use pointers for flexibility */
    template <typename T>
    class LoggerDistributor : public ChainLogger<T>
    {
    public:
        typedef typename ChainLogger<T>::SampleType SampleType;
        typedef typename ChainLogger<T>::Generator  Generator;
        typedef typename ChainLogger<T>::Evaluator  Evaluator;
        typedef typename ChainLogger<T>::Logger     Logger;

    public:
        typedef std::vector<Logger*> LogList;
        typedef typename LogList::iterator LogIterator;

    public:
        LoggerDistributor() :
                Destinations()
        {}

        /// dtor does nothing, clean up memory yourself
        virtual ~LoggerDistributor() {}

        virtual void notifyAccept( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
        {
            for( LogIterator it = Destinations.begin(); it != Destinations.end(); it++ )
                (*it)->notifyAccept( cparm, pValue, proposal, evaluator );
        }
        virtual void notifyReject( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
        {
            for( LogIterator it = Destinations.begin(); it != Destinations.end(); it++ )
                (*it)->notifyReject( cparm, pValue, proposal, evaluator );
        }
        virtual void notifyReset( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
        {
            for( LogIterator it = Destinations.begin(); it != Destinations.end(); it++ )
                (*it)->notifyReset( cparm, pValue, proposal, evaluator );
        }

        void add ( Logger* newlogger )
        {
            Destinations.push_back( newlogger );
        }

    protected:
        LogList Destinations;
    };

}

#endif   /* ----- #ifndef LOGGERDISTRIBUTOR_INC  ----- */
