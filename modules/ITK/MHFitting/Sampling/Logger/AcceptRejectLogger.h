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

#ifndef  __ACCEPT_REJECT_LOGGER
#define  __ACCEPT_REJECT_LOGGER

#include <vector>

#include "../MarkovChain.h"

namespace sampling
{
    /** \brief Copy logging stream to multiple loggers - Use pointers for flexibility */
    template <typename T>
    class AcceptRejectLogger : public ChainLogger<T>
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
        AcceptRejectLogger(std::ofstream& logstream) :
                m_logstream(logstream), m_numAccepted(0), m_numRejected(0)
        {}


        virtual ~AcceptRejectLogger() {   }

        virtual void notifyAccept( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
        {
            m_numAccepted += 1;
            m_logstream << "accept" << "," << proposal->getName() << "," << evaluator->evalSample(cparm) << "," << pValue << "," << m_numAccepted / static_cast<double>(m_numAccepted + m_numRejected)  << std::endl;
        }
        virtual void notifyReject( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
        {
            m_numRejected += 1;
            m_logstream << "reject" << "," << proposal->getName() << "," << evaluator->evalSample(cparm) << "," << pValue << "," << m_numAccepted / static_cast<double>(m_numAccepted + m_numRejected) << std::endl;
        }
        virtual void notifyReset( const T& cparm, const double& pValue, Generator* proposal,  Evaluator* evaluator )
        {
        }


    protected:
        std::ofstream& m_logstream;
        unsigned m_numAccepted;
        unsigned m_numRejected;
    };

}

#endif
