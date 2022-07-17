#ifndef __SUPERTILES_PLACE_LOG__
#define __SUPERTILES_PLACE_LOG__

#include "supertiles_place_term.h"
#include "supertiles_place_plan.h"
#include <sstream>

namespace supertiles
{
  namespace place
  {
    struct PassTimings
    {
      friend std::ostream& operator<< (std::ostream &out, const PassTimings &e)
      {
	out << e.total << ","<<e.plan << ","<<e.exchangeAdapt;
	return out;
		  
      }
      double total=0.;
      double plan=0.;
      double exchangeAdapt=0.;
    };
    
    template<typename PASS_TIMINGS>
    struct LogElem
    {
      using D=double;
      using I=int64_t;

      I passN=-1;
      D timestamp=-1.;
      D cost=-1.;
      uint32_t change=false;
      I plan_level=-1;
      bool plan_aggregateExchangeMode=false;
      PASS_TIMINGS pt;
      
      friend std::ostream& operator<< (std::ostream &out, const LogElem &e)
      {
	out
	  << e.passN << ","
	  << e.timestamp << ","
	  << e.cost << ","
	  << e.change << ","
	  << e.plan_level << ","
	  << e.plan_aggregateExchangeMode << ","
	  << e.pt;
	return out;
      }
      
    };
    
    template<typename TERM,typename PASS_TIMINGS>
    struct Log
    {
      using D=double;
      using I=int64_t;
      
      Log(TERM& term_in) :
	term(term_in)
      {}	
      
      void operator()(D cost, uint32_t change, const uint32_t level, const bool aggregateExchangeMode, const PASS_TIMINGS pt)
      {
	LogElem<PASS_TIMINGS> e;

	e.passN = term.passN;
	e.timestamp = term.timer.get_s();

	e.cost=cost;
	e.change=change;

	e.plan_level=level;
	e.plan_aggregateExchangeMode=aggregateExchangeMode;

	e.pt=pt;
	l.push_back(e);
      }

      void writeCSV(const std::string fname, bool append) const
      {
	std::cout << "write log to " << fname << std::endl;
	std::stringstream ss;
	for(const auto & e : l)
	  ss << e << std::endl;
	std::ofstream fout(fname.c_str(), append ? std::ios_base::app : std::ios_base::trunc);
	fout << ss.str();
	fout.close();
      }

      const TERM& term;
      std::vector<LogElem<PASS_TIMINGS>> l;      
    };
    
  }
}
#endif //__SUPERTILES_PLACE_LOG__
