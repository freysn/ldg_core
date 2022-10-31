#ifndef __SUPERTILES_PLACE_TERM__
#define __SUPERTILES_PLACE_TERM__

#include "helper_signal.h"
#include "helper_ChronoTimer.h"

namespace supertiles
{
  namespace place
  {
    template<typename D>
    struct Term
    {
      Term(/*bool useSignalHandler=true*/)
      {
	//if(useSignalHandler)
	  initSignalHandler();
      }

      template<typename PO>
      Term (const PO& po)
      {
	maxNSeconds=po.termTime;
	maxNPasses=po.termIterations;
	maxNSame=po.termSame;
	initSignalHandler();
      }

      void initSignalHandler()
      {
	signalHandler = new helper::SignalHandler();
      }

      ~Term()
      {
	// if(signalHandler)
	//   delete signalHandler;
      }
      
      bool operator()(bool change, bool justChecking=false)
      {
	if(!justChecking)
	  passN++;
	
	if(signalHandler && signalHandler->gotExitSignal())
	  return true;
	
	bool term = (passN>=maxNPasses);

	if(term)
	  return true;
    
	
	if(!justChecking)
	  {
	    if(change)
	      sameCnt=0;
	    else
	      sameCnt++;
	  }
	if(sameCnt >= maxNSame)
	  return true;

	//std::cout << "timer.get_s() " << timer.get_s() << " vs " << maxNSeconds << std::endl;
	if(timer.get_s() >= maxNSeconds)
	    return true;
	    
    
	return false;
      }

      auto getPassN() const
      {
	return passN;
      }

      size_t passN=0;

      size_t sameCnt=0;

      size_t maxNPasses=std::numeric_limits<size_t>::max();
      size_t maxNSame = std::numeric_limits<size_t>::max();;

      double maxNSeconds = std::numeric_limits<double>::max();

      helper::SignalHandler* signalHandler=0;
      helper::ChronoTimer timer;

      double cost=-1;
    };    
  }
}
#endif //__SUPERTILES_PLACE_TERM__
