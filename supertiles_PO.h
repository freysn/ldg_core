#ifndef __SUPERTILES_PO__
#define __SUPERTILES_PO__

#include <vector>
#include "helper_volData/splitStr.h"
#include "helper_lexicalCast.h"
#include "supertiles_configTypes.h"
#include <iostream>
#include <cassert>
#include "helper_assert.h"

  enum dataInputFormat_t
  {
    dif_spatial=0,
    dif_feature=1
  };

  enum distMode_t
  {
    dm_euclidean2=0,
    dm_euclidean=1
  };


template<typename E>
E getEnumType(const std::string s)
{
  return static_cast<E>(std::stoi(s));
}
  
  class supertiles_PO
  {
  public:
    supertiles_PO() :
    termTime(100000000.),
      termImprovement(0.000001),      
      termIterations(std::numeric_limits<size_t>::max()),
    
    termSame(std::numeric_limits<size_t>::max()),
    distanceLimit(0.000001),
      preNormData(false),
    
    //outDir("out_interp/"),
    outDir(""),
      //initNearest(true),
      alignCenters(false),
    nIntermediates(10),
      doWriteResult(true),
      doWriteLog(true),
      maxValue(255),
      dataInputFormat(dif_spatial),
    distMode(dm_euclidean2)
    //,
    //binarize(-1.)
	{	  
	};
    /*
    static constexpr unsigned int str2int(const char* str, int h = 0)
    {
      return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
    }
    */

    template<typename ARGV>
    bool handleProgramOptions(int argc, const ARGV argv, bool doHandleVM=true)
    {      
      
      for(int argi=1; argi < argc; argi++)
	{
	  
	  auto nextStr = [&argi, argc, argv]()
	    {
	      argi++;
	      hassert(argi < argc);
	      return std::string(argv[argi]);
	    };

	  std::string __s = std::string(argv[argi]);

	  auto trim = [](std::string s)
	    {
	      if(s=="")
		return s;
	      
	      const auto begin = s.find_first_not_of(" ");
	      if(begin == std::string::npos)
		return std::string("");
	      
	      const auto end = s.find_first_of(" ", begin);
	      size_t len = std::string::npos;
	      if(end != std::string::npos)
		len = end-begin;
	      //std::cout << "|" << s << "|" << begin << "|" << end << "|" << len << "|\n";
	      return s.substr(begin, len);
	    };

	  auto c = [__s, trim](std::string x, std::string m="")
	    {

	      // print help message
	      if(m != "" && (__s == "-h" || __s == "--help"))
		  std::cout << "help: " << x << " | " << m << std::endl;
	      
	      std::vector<std::string> xv = split(x, ',');	      
	      
	      for(auto e : xv)
		{
		  const bool match = (__s == trim(e));

		  //std::cout << __s << " vs " << trim(e) << " " << match << std::endl;
		  if(match)
		    return true;
		}
	      return false;
	    };

	  if(c("-d, --datconfs", "data config file")) {configFNames.push_back(nextStr());}
	  else if(c("-r, --rep"))
	    {repFNames.push_back(nextStr());}
	  else if(c("-t, --trans", "transfer function file")){transFuncFNames.push_back(nextStr());}
	  else if(c("-s, --steps")){timeSteps.push_back(helper::lexicalCast<int>(nextStr()));}
	  else if(c(", --maxValue")){helper::lexicalCast(maxValue, nextStr());}
	  
	  else if(c("--distFuncType"))
	    {distFuncType = getEnumType<distFuncType_t>(nextStr());}
	  // else if(c(", --filterType"))
	  //   {filterType = getEnumType<filterType_t>(nextStr());}
	  else if(c("--distVisType"))
	    {distVisType = getEnumType<distVisType_t>(nextStr());}
	  else if(c("--repAggregationType"))
	    {repAggregationType = getEnumType<repAggregationType_t>(nextStr());}
	  
	  else if(c("--nIntermediates", "number of output intermediates")){helper::lexicalCast(nIntermediates, nextStr());}
	  else if(c("--termTime"))
	    {
	      helper::lexicalCast(termTime, nextStr());
	    }
	  else if(c("--nNeighbors"))
	    {
	      helper::lexicalCast(nNeighbors, nextStr());
	    }
	  else if(c("--neighborFac"))
	    {
	      helper::lexicalCast(neighborFac, nextStr());
	    }
	  else if(c("--planChunkSize"))
	    {
	      helper::lexicalCast(planChunkSize, nextStr());
	    }
	  else if(c("--termImprovement"))
	    {
	      helper::lexicalCast(termImprovement, nextStr());
	    }
	  else if(c("--maxNLevelsUp"))
	    {
	      helper::lexicalCast(maxNLevelsUp, nextStr());
	    }
	  else if(c("--maxNodeExchangeLevel"))
	    {
	      helper::lexicalCast(maxNodeExchangeLevel, nextStr());
	    }
	  else if(c("--termIterations"))
	    {
	      helper::lexicalCast(termIterations, nextStr());
	    }
	  else if(c("--termSame"))
	    {
	      helper::lexicalCast(termSame, nextStr());
	    }
	  else if(c("--distancelimit"))
	    {
	      helper::lexicalCast(distanceLimit, nextStr());
	    }
	  else if(c("--batchSize"))
	    {
	      helper::lexicalCast(batchSize, nextStr());
	    }
	  else if(c("--seed"))
	    {
	      helper::lexicalCast(seed, nextStr());
	    }
	  else if(c("--nTilesAssign"))
	    {
	      helper::lexicalCast(nTilesAssign, nextStr());
	    }
	  else if(c("--preNormData"))
	    {preNormData=true;}
	  // else if(c("-r,--reductionFactor", "reductionFactor"))
	  //   {reductionFactors.push_back(helper::lexicalCast<int>(nextStr()));}
	  
	  else if(c("--selectSliceParams", "selectSlice"))
	    {selectSliceParams.push_back(helper::lexicalCast<double>(nextStr()));}
	  // else if(c("--binarize", "selectSlice"))
	  //   {binarize = helper::lexicalCast<double>(nextStr());}
	  // else if(c("--blurRadiusLeaf"))
	  //   {blurRadiusLeaf = helper::lexicalCast<double>(nextStr());}
	  else if(c("--outputImgScale"))
	    {outputImgScale = helper::lexicalCast<double>(nextStr());}
	  else if(c("--outDir", "output directory"))
	    {
	      outDir = nextStr();	      
	    }
	  else if(c("--writeIntermediateAssignments"))
	    {
	      writeIntermediateAssignments = nextStr();
	    }
	  else if(c("--writeGridAssignmentTXT"))
	    {
	      writeGridAssignmentTXT = nextStr();
	    }
	  else if(c("--alignCenters", "align centers"))
	    {alignCenters = true;}
	  else if(c("--noImgOut"))
	    {imgOut = false;}
	  else if(c("--imgOutOpts"))
	    {imgOutOpts = nextStr();}
	  else if(c("--noDisparityIndicator"))
	    {drawDisparityIndicator = false;}
	  else if(c("--imgOut"))
	    {imgOut = true;}
	  else if(c("--getCost"))
	    {getCost = true;}
	  else if(c("--doNotWriteResult", "do not write result"))
	    {doWriteResult = false;}
	  else if(c("--runAgents", "runAgents"))
	    {runAgents = nextStr();}
	  else if(c("--doNotWriteLog", "do not write log"))
	    {doWriteLog = false;}
	  else if(c("--drawLevels"))
	    {drawLevels = true;}
	  else if(c("--drawInputTiles"))
	    {drawInputTiles = true;}
	  else if(c("--termDistRef"))
	    {termDistRef = nextStr();}
	  else if(c("--repCacheFName"))
	    {repCacheFName = nextStr();}
	  else if(c("--writeDisparitesOnly"))
	    {writeDisparitesOnly = nextStr();}
	  else if(c("--rebuildRepCache"))
	    {rebuildRepCache = true;}
	  else if(c("--dataInputFormat"))
	    {dataInputFormat = static_cast<dataInputFormat_t>
		(helper::lexicalCast<int>(nextStr()));}
	  else if(c("--dataInputFormat"))
	    {dataInputFormat = static_cast<dataInputFormat_t>
		(helper::lexicalCast<int>(nextStr()));}
	  else if(c("-a, --loadAssignments"))
	    {loadAssignments = nextStr();}
	  else if(c("--timeSteps"))
	    {timeSteps.push_back(helper::lexicalCast<size_t>(nextStr()));}
	  else if(c("-h, --help", "print help message")) {}
	  else
	    {
	      std::cerr << __PRETTY_FUNCTION__
			<< " argument not recognized: " << argv[argi] << std::endl;
	      assert(false);
	      exit(-1);
	    }
	}

      if(timeSteps.empty())
	timeSteps.resize(1,0);

      if(doHandleVM)
	handleVM();

      return true;
    }
        

    void handleVM()
    {
      if (configFNames.empty())
	{
	  // configFNames.push_back("resources/i0.config");
	  // configFNames.push_back("resources/i1.config");
	  configFNames.push_back("65536");
	}

      if(timeSteps.empty())
	{
	  timeSteps.push_back(0);
	}
    }

    //std::vector<int> reductionFactors;
    std::vector<double> selectSliceParams;
    std::vector<std::string> configFNames;
    std::vector<std::string> repFNames;
    std::vector<std::string> transFuncFNames;
    std::vector<size_t> timeSteps;
    double termTime;
    double termImprovement;
    size_t termIterations;
    size_t termSame;
    double distanceLimit;
    bool preNormData;
    std::string outDir;
    //bool initNearest;
    bool alignCenters;

    size_t nIntermediates;

    std::string termDistRef;

    bool doWriteResult;
    bool doWriteLog;
    size_t maxValue;

    dataInputFormat_t dataInputFormat;
    distMode_t distMode;

    //double binarize;

    //double blurRadiusLeaf = 2.;

    std::string loadAssignments="";

    bool drawInputTiles=false;
    bool drawLevels=false;

    distFuncType_t distFuncType=distFuncType_none;
    //filterType_t filterType=filterType_box;
    distVisType_t distVisType=distVisType_none;
    repAggregationType_t repAggregationType=repAggregationType_average;

    bool drawScoresSelection=true;

    int nNeighbors=4;
    //double neighborFac=1./8.;
    double neighborFac=1./4.;

    size_t planChunkSize=6;
    size_t maxNLevelsUp=-1;
    size_t maxNodeExchangeLevel=-1;

    std::string writeIntermediateAssignments;

    std::string writeGridAssignmentTXT;    

    size_t seed=1337;
    size_t batchSize=128;

    double outputImgScale=.1;
    bool imgOut=true;

    bool drawDisparityIndicator=true;

    size_t nTilesAssign=-1;

    std::string imgOutOpts;

    std::string repCacheFName;
    bool rebuildRepCache=false;

    std::string writeDisparitesOnly;

    bool getCost=false;

    std::string runAgents;
  };

template<typename ARGV>
auto initPO(int argc, const ARGV argv)
{
  supertiles_PO po;
  if(!po.handleProgramOptions(argc, argv, false))
    {
      std::cout << "cannot handle program options ... exiting\n";
      exit(-1);
    }
  return po;
}


#endif //# __SUPERTILES_PO__
