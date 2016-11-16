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
 *       Filename:  RandomGenerator.h
 *
 *    Description:  A class for Random Number Generation, based on Boost
 *
 *        Version:  1.0
 *        Created:  20.04.2010 09:03:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */


#ifndef  RANDOMGENERATOR_INC
#define  RANDOMGENERATOR_INC

#include <vector>

#include <boost/random.hpp>

namespace sampling
{

  //TODO not really a lib, a bunch of header files
  //TODO different distributions, fast access
  //TODO iterator access, stream access!
  //

  class NormalFunctor;

  /** \brief Functor interface */
  class RandomFunctor
  {
    public:
      virtual double operator() () = 0;
  };

  /** \brief A Random Generator class, as a wrapper around boost::mt19937
   *
   *  It can be used to solve the seeding problem and provide easier access to distributions
   *  Direct conversions between instances of this class and GeneratorType are possible */
  class RandomGenerator
  {
    public:
      typedef boost::mt19937 GeneratorType;

      typedef RandomGenerator* Ptr;

    protected:
      GeneratorType generator;

      boost::uniform_real<double> distUniformDbl;
      boost::normal_distribution<double> distNormalDbl;

      boost::variate_generator<GeneratorType&, boost::normal_distribution<double> > genNormal;

    public:
      /** \brief Construct from GeneratorType */
      RandomGenerator( const GeneratorType& rndEngine ) : generator( rndEngine ), distUniformDbl(0.,1.), distNormalDbl(), genNormal( generator, distNormalDbl ) {}
      /** \brief Construct with explicit seed */
      RandomGenerator( unsigned int seed ) : generator( seed ), distUniformDbl(0.,1.), distNormalDbl(), genNormal( generator, distNormalDbl ) {}
      /** \brief Construct with seed from getGoodRandomSeed */
      RandomGenerator() : generator( RandomGenerator::getGoodRandomSeed() ), distUniformDbl(0.,1.), distNormalDbl(), genNormal( generator, distNormalDbl ) {}

      /** \brief Get a uniform double in [0,1) */
      double uniformDbl()
      {
        return distUniformDbl(generator);
      }

      /** \brief Get a uniform integer */
      double uniformDbl( double lowerBound, double upperBound )
      {
        return uniformDbl()*(upperBound-lowerBound)+lowerBound;
      }

      /** \brief Get a standard normal double */
      double normalDbl()
      {
        return genNormal();
      }

      /** \brief Get a uniform integer */
      int uniformInt( int lowerBound, int upperBound )
      {
        return int( uniformDbl()*(upperBound-lowerBound)+lowerBound );
      }

      std::vector<int> uniformInts( int N, int lowerBound, int upperBound );

      GeneratorType& getEngine()
      {
        return generator;
      }

      /** \brief Provides direct conversion to GeneratorType */
      operator GeneratorType& ()
      {
        return getEngine();
      }
      /** \brief Provides direct conversion to GeneratorType */
      operator GeneratorType()
      {
        return GeneratorType(getEngine());
      }

      /** \brief assignment: std::copy state of generator */
      void operator= ( const RandomGenerator& that )
      {
        generator = GeneratorType( that.generator );
      }

      size_t operator()( size_t upperBound )
      {
        return uniformInt( 0, upperBound );
      }

      double operator()()
      {
        return uniformDbl();
      }

      /** \brief Generate templated numeric type numbers */
      template <typename nT>
      std::vector<nT> normals( int N, double mean, double variance )
      {
        std::vector<nT> rndPos;
        rndPos.reserve(N);
        for ( int i = 0; i < N; i++ )
        {
          rndPos.push_back( static_cast<nT>( normalDbl()*variance + mean ) );
        }
        return rndPos;
      }

      /** \brief Generate a random permutation of N elements */
      std::vector<size_t> randomPermutation( const size_t n );

    private:

      // static
    public:
      /** \brief Returns a random seed (unsigned int) from /dev/urandom - linux dependent! */
      static unsigned int getGoodRandomSeed();
  };


  // TODO Functors: inherit from RandomGenerator or leave it as it is?
  // TODO what about boost::variate_generator< GeneratorType&, DISTRIBUTION > ??

  /** \brief Generate random normal doubles - overloaded ()-operator available */
  class NormalFunctor : RandomFunctor
  {
    public:
      NormalFunctor( RandomGenerator& rndEngine ) : generator( rndEngine ) {}

      double operator() (double mean, double stddev )
      {
        return generator.normalDbl()*stddev + mean;
      }
      double operator() ()
      {
        return generator.normalDbl();
      }

      // read a double
      void operator >> ( double& target )
      {
        target = operator()();
      }
    private:
      RandomGenerator& generator;
  };

  /** \brief Generate random uniform doubles - overloaded ()-operator available */
  class UniformFunctor : RandomFunctor
  {
    public:
      UniformFunctor( RandomGenerator& rndEngine ) : generator( rndEngine ) {}

      double operator() (double lower, double upper )
      {
        return generator.uniformDbl()*(upper-lower)+lower;
      }
      double operator() ()
      {
        return generator.uniformDbl();
      }

      // read a double
      void operator >> ( double& target )
      {
        target = operator()();
      }
    private:
      RandomGenerator& generator;
  };

}
#endif   /* ----- #ifndef RANDOMGENERATOR_INC  ----- */
