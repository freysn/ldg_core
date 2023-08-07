#ifndef __SUPERTILES_PLACE_DRAWOPTS__
#define __SUPERTILES_PLACE_DRAWOPTS__

#include "helper_color.h"

namespace supertiles
{
  struct DrawOpts
  {
    DrawOpts(const std::string& imgOutOpts)
    {
      parse(imgOutOpts);
    }

    void parse(const std::string& imgOutOpts)
    {
      const auto keyValuePairs = helper::split(imgOutOpts,',');
      for(const auto & e : keyValuePairs)
	{
	  if(e.empty())
	    continue;
	  const auto keyValue=helper::split(e,':');
	  hassertm(keyValue.size()==2, e);
	  if(keyValue[0]=="import")
	    {
	      std::vector<std::string> lines;
	      const bool success=
		helper::readFileAuto(lines, keyValue[1], "txt");

	      std::cout << "successfully loaded imported file "
			<< keyValue[1] << ": " << success
			<< std::endl;
	      
	      for(const auto & e : lines)
		{
		  assert(e != imgOutOpts);
		  parse(e);
		}
	    }
	  else if(keyValue[0]=="level" || keyValue[0]=="height")
	    levels.push_back(std::stoi(keyValue[1]));	  
	  else if(keyValue[0]=="disparity")
	    disparities.push_back(std::stod(keyValue[1]));
	  else if(keyValue[0]=="annotateTiles")
	    annotateTiles=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="pdf")
	    do_pdf=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="png")
	    do_png=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="sizeLegend")
	    do_drawSizeLegend=true;
	  else if(keyValue[0]=="maxNTiles")
	    maxNTiles=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="swapXY")
	    swapXY=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="level0border")
	    level0Border=std::stod(keyValue[1]);
	  else if(keyValue[0]=="levelThresholdMin")
	    levelThresholdMin=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="zoomActiveNodes")
	    zoomActiveNodes=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="fullVoidRatioAdapt")
	    fullVoidRatioAdapt=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="shownNodes_borderCol")
	    shownNodes_borderCol=helper::hex2rgba(keyValue[1]);
	  else if(keyValue[0]=="shownNodes_ann_borderCol")
	    shownNodesv_ann_borderCol.push_back(helper::hex2rgba(keyValue[1]));
	  else if(keyValue[0]=="background")
	    background=helper::hex2rgba(keyValue[1]);
	  else if(keyValue[0]=="disparityIndicatorCol")
	    disparityIndicatorCol=helper::hex2rgba(keyValue[1]);
	  else if(keyValue[0]=="blur")
	    blur=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="mcmc_cm")
	    mcmc_cm=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="do_drawCM")
	    do_drawCM=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="do_uniformBorder")
	    do_uniformBorder=std::stoi(keyValue[1]);
	  else if(keyValue[0]=="tileRepDimScale")
	    tileRepDimScale=std::stod(keyValue[1]);
	  else if(keyValue[0]=="shownNodes")
	    {
	      //std::vector<std::string> fnames
	      //=helper::genFileNames(std::vector<std::string>(1,keyValue[1]));

	      std::vector<std::string> fnames = helper::cmd_ls(keyValue[1]);
	      std::sort(fnames.begin(), fnames.end());
	      shownNodes.resize(fnames.size());
	      int cnt=0;
	      for(const auto & e : fnames)
		{
		  std::cout << "load shownNodes " << e << std::endl;
		  const bool success
		    =helper::readFileAuto(shownNodes[cnt++], e);
		  hassertm(success, /*keyValue[1]*/e);
		}
	      // hassertm2(shownNodes.size()==qt.nElems(),
	      // 		shownNodes.size(),
	      // 		qt.nElems());
	    }
	  else if(keyValue[0]=="shownNodes_ann")
	    {
	      std::vector<uint8_t> shownNodes_ann;
	      const bool success
		=helper::readFileAuto(shownNodes_ann, keyValue[1]);
	      hassert(success);
	      shownNodesv_ann.push_back(shownNodes_ann);
	      // hassertm2(shownNodes_ann.size()==qt.nElems(),
	      // 		shownNodes_ann.size(),
	      // 		qt.nElems());
	    }
	  else if(keyValue[0]=="borderLineWidthScale")
	    borderLineWidthScale=std::stod(keyValue[1]);
	  else if(keyValue[0]=="levelIndicatorScale")
	    levelIndicatorScale=std::stod(keyValue[1]);
	  else if(keyValue[0]=="disparityIndicatorWidthScale")	    
	    disparityIndicatorWidthScale=std::stod(keyValue[1]);
	  else if(keyValue[0]=="level0BorderScale")	    
	    level0BorderScale=std::stod(keyValue[1]);
	  else if(keyValue[0]=="individualTileImgs")
	    individualTileImgs=keyValue[1];
	  else if(keyValue[0]=="chartLineWidthScale")
	    chartLineWidthScale=std::stod(keyValue[1]);
	  else if(keyValue[0]=="shownNodes_inactiveOpMod")
	    shownNodes_inactiveOpMod=std::stod(keyValue[1]);
	  else if(keyValue[0]=="mcmc_labelsFName")
	    mcmc_labelsFName=keyValue[1];
	  else if(keyValue[0]=="shownNodes_ghostArea")
	    {
	      const auto valuePair=helper::split(keyValue[1],'|');
	      assert(valuePair.size()==2);
	      if(valuePair.size()==2)		
		shownNodes_ghostArea=V2<double>(std::stod(valuePair[0]), std::stod(valuePair[1]));
	    }	  	  
	  else
	    {
	      std::cerr << "Error: keyword " << keyValue[0]
			<< " not recognized, input string: |" <<  imgOutOpts << "|\n";
	      throw "something";
	    }
	}
    }

    std::vector<size_t> levels;
    std::vector<double> disparities;
        
    bool do_pdf=true;
    bool do_png=true;
    size_t maxNTiles=128*1024;
    bool swapXY=false;

    double level0Border=-1.;
    uint32_t levelThresholdMin=-1;

    //----------
    std::vector<std::vector<uint8_t>> shownNodes;
    V4<double> shownNodes_borderCol = V4<double>(0, 0,0, 0);
    std::vector<std::vector<uint8_t>> shownNodesv_ann;
    std::vector<V4<double>> shownNodesv_ann_borderCol;
    //= V4<double>(0, 0, 0, 0);
    bool zoomActiveNodes=false;
    //V4<double> disparityIndicatorCol=V4<double>(.3, .3, .3, 1.);
    V4<double> disparityIndicatorCol=V4<double>(0., 0., 0., 0.);
    bool fullVoidRatioAdapt=true;
    double borderLineWidthScale=1.;
    double levelIndicatorScale=1.;
    double disparityIndicatorWidthScale=3.;
    double level0BorderScale=1.;
    V4<double> background=V4<double>(0., 0., 0., 0.);
    bool blur=false;
    std::string individualTileImgs;

    double shownNodes_inactiveOpMod=.3;
    V2<double> shownNodes_ghostArea=V2<double>(.15, .15);

    double chartLineWidthScale=1.;

    double tileRepDimScale=1.;
    bool do_drawSizeLegend=false;
    
    uint32_t mcmc_cm=0;
    bool do_drawCM=true;

    bool annotateTiles=false;

    bool do_uniformBorder=false;

    std::string mcmc_labelsFName;
  };
}

#endif // __SUPERTILES_PLACE_DRAWOPTS__
