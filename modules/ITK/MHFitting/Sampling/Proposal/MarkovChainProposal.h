/*
 * =====================================================================================
 *
 *       Filename:  MarkovChainProposal.h
 *
 *    Description:  Proposal from a Markov Chain (can be used to stack them hierarchically)
 *
 *        Version:  1.0
 *        Created:  01/14/2011 09:40:23 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sandro Schönborn (ses), sandro.schoenborn@unibas.ch
 *        Company:  University of Basel
 *
 * =====================================================================================
 */

#ifndef  MARKOVCHAINPROPOSAL_INC
#define  MARKOVCHAINPROPOSAL_INC

#include "../MarkovChain.h"

namespace sampling
{

  template <typename T>
  class MarkovChainProposal : public ProposalGenerator<T>
  {
    public:
      typedef MarkovChain<T> MC;
      typedef DistributionEvaluator<T> Evaluator;

    public:
      MarkovChainProposal( MC* pChain, size_t nsubsteps, bool reset = false ) :
        pMC( pChain ),
        pEval( 0 ),
        nSteps( nsubsteps ),
        bReset( reset )
      {
      }

      MarkovChainProposal( MC* pChain, Evaluator* transition_evaluator, size_t nsubsteps, bool reset = false ) :
        pMC( pChain ),
        pEval( transition_evaluator ),
        nSteps( nsubsteps ),
        bReset( reset )
      {
      }

      void generateProposal( T& rSample, const T& rCurrent )
      {
        if ( bReset )
          pMC->setState( rCurrent ); // restart chain at rCurrent

        pMC->current( rSample );
        for ( size_t i = 0; i < nSteps; i++ )
          pMC->next( rSample );
      }

      double transitionProbability( const T& from, const T& to )
      {
        // hmmmmmm ?????
        // Assume a working chain as base -> Q(from,to) = Q(to) = P_chain(to)
        // Ratio is built by corresponding algorithm
        if ( pEval )
          return pEval->evalSample( to );
        else
          return 0.0;
      }

    private:
      MC* pMC;
      Evaluator* pEval;
      size_t nSteps;
      bool bReset;
  };

}
#endif   /* ----- #ifndef MARKOVCHAINPROPOSAL_INC  ----- */


