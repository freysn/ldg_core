#ifndef __SUPERTILES_PLACE_DATA__
#define __SUPERTILES_PLACE_DATA__

#define NO_CUDA_TYPES
#include "helper_volData/UserData.h"

#include "helper_cmd.h"
#include "helper_readDataFromConfigs.h"
#include "helper_imath.h"
#include "helper_io.h"

namespace supertiles
{
  namespace place
  {

    template<typename D>
    D genSynthDataElem(const V2<size_t> /*gridId*/, const V2<size_t> /*gridDim*/,
		       const V2<size_t> /*tileElemId*/, const V2<size_t> /*tileDim*/)
    {
      std::cerr << "THIS SHOULD NEVER BE INSTANTIATED\n";
      assert(false);
      exit(-1);
    }

    template<>
    double genSynthDataElem
    <double>
    (const V2<size_t> gridId, const V2<size_t> gridDim,
     const V2<size_t> /*tileElemId*/, const V2<size_t> /*tileDim*/)
    {
      return std::max(static_cast<double>(gridId.x)/gridDim.x,
		      static_cast<double>(gridId.y)/gridDim.y);
    }

    template<>
    V4<unsigned char> genSynthDataElem
    <V4<unsigned char>>
    (const V2<size_t> gridId, const V2<size_t> gridDim,
     const V2<size_t> /*tileElemId*/, const V2<size_t> /*tileDim*/)
    {
      return V4<unsigned char>(gridId.x*255./gridDim.x,
			       48,
			       gridId.y*255./gridDim.y,
			       255
			       );
    }

    auto readUserData(std::string configFName)
    {
      using UserData_t = UserData<V3<int>,V3<double>>;
      UserData_t ud;
      {
	const bool success = ud.readConfigOnly(configFName.c_str());
	hassert(success);
      }
      return ud;
    }

    // template<typename F>
    // auto loadDataConfig(const std::string configFName)
    // {

    // }

    template<typename F>
    auto loadData(const std::string basePath,
		  const size_t firstNFiles=std::numeric_limits<size_t>::max(),
		  bool normalize=true,
		  bool force2d=false)
    {
      std::vector<F> data;

      V4<size_t> dim(0,0,0,0);

      //if(basePath != "")
      assert(basePath != "");
	{
	  std::cout << "load config files with pattern |"
		    << basePath << "|\n";

	  std::vector<std::string> configFiles;

#ifndef NO_CMD
	  configFiles= helper::cmd_ls(basePath);
#else
	    configFiles=std::vector<std::string>(1, basePath);
#endif
	    //std::cout << "there are << " << configFiles.size() << " config files\n";
	  std::sort(configFiles.begin(), configFiles.end());
	  configFiles.resize(std::min(firstNFiles, configFiles.size()));

	  size_t nMembers = 0;
	  for(const auto & configFile : configFiles)
	    {
	      std::cout << "load file " << configFile << std::endl;
	      const auto origData = helper::readDataFromConfig<F>(configFile);
	      dim = std::get<1>(origData);
	      assert(!std::get<0>(origData).empty());
	      //for(const auto timeStep : std::get<0>(origData))
	      size_t cnt=0;
	      for(size_t t=0; t<std::get<0>(origData).size(); t++)
		{
		  auto timeStep = std::get<0>(origData)[t];
		  data.insert(data.end(), timeStep.begin(), timeStep.end());
		  cnt++;
		}
	      nMembers += cnt;
	    }
	  //dim.w = std::min(nMembers, firstNFiles);
	  dim.w = nMembers;

	  std::cout << "dim: " << dim << std::endl;

	  {
	    const auto nElemsDataConfig=helper::iiii2n(dim);
	    if(nElemsDataConfig>data.size())
	      {
		std::cerr << "ERROR: expected " << nElemsDataConfig << " elements, but input data only features " << data.size() << std::endl;
		throw "size of provided data smaller than defined in config file";
	      }
	    else if(nElemsDataConfig < data.size())
	      {
		std::cout << "data size is " << data.size() << " but config file states " << nElemsDataConfig << " elems, resizing data\n";
		data.resize(nElemsDataConfig);
	      }
	  }


	  if(normalize)
	    {
	      std::cout << "normalize data\n" << std::endl;
	      helper::normalize_minmax(data.begin(), data.end());
	    }

	  if(force2d)
	    {
	      dim.w *= dim.z;
	      dim.z = 1;
	    }
	}
#if 0
      else
	{
	  std::cout << "generate test data\n";
	  dim.w = 8;
	  dim.z=1;

	  dim.x=dim.y=16;

	  const size_t n = helper::iiii2n(dim);

	  data.resize(n);

	  std::mt19937 gen(1337);

	  std::uniform_real_distribution<> dis(0.0, 1.0);

	  for(auto & e : data)
	    e = dis(gen);

	  auto art = [&](auto off)
	  {
	    for(const auto m : {2,4})
	      for(const auto y : helper::range_bn(off.y, 8))
		for(const auto x : helper::range_bn(off.x, 8))
		  data[helper::iiii2i(V4<size_t>(x,y,0,m), dim)] = 1.;
	  };

	  art(V2<size_t>(0,0));
	  art(V2<size_t>(8,8));
	}
#endif
      return std::make_tuple(data, dim);
    }

    template<typename DIM>
    auto nTilesFromConfigDim(const DIM& dim)
    {
      return dim.z;
    }

    template<typename DIM>
    auto initGridNTiles(size_t nFull, size_t nTilesAssign)
    {
      DIM gridDim;
      size_t nTiles=nFull;
      if(nTilesAssign!=static_cast<size_t>(-1))
	{
	  assert(nFull<=nTilesAssign);
	  nFull=nTilesAssign;
	}

      const auto exp = helper::ilog4_ceil(nFull);
      gridDim.x = helper::ipow2(exp);
      assert(gridDim.x==std::pow(2, exp));
      nTiles=nFull;
      gridDim.y = gridDim.x;


      return std::make_tuple(gridDim, nTiles);
    }

    template<typename D, typename PO>
    auto initData(const PO& po, std::string configFName)
    {
      V2<size_t> gridDim(16, 16);

      // use loaded image to initalize data values and compute distances

      std::vector<D> tileData;
      V2<int> tileDim(1,1);

      size_t nTiles=-1;
      if(!po.configFNames.empty())
	{
	  assert(po.configFNames.size()==1);
	  const std::string fname(configFName);
	  std::cout << "read data from config file " << configFName << std::endl;

	  //enum class dataMode_t {signal, image};

	  //dataMode_t dm = dataMode_t::image;

	  // if(fname.find("5Jets") != std::string::npos)
	  // 	dm = dataMode_t::signal;

	  //const bool doNormalize = (dm != dataMode_t::signal);

	  auto [data, dim] = loadData<D>(configFName,
					 std::numeric_limits<size_t>::max(),
					 po.preNormData);

	  size_t nFull=0;

	  {
	    //dm = dataMode_t::image;

	    //dm = dataMode_t::signal;
	    tileDim.x=dim.x;
	    tileDim.y=dim.y;

	    //nFull = dim.z;
	    nFull=nTilesFromConfigDim(dim);
	    assert(dim.w==1);
	  }


	  std::tie(gridDim, nTiles) = initGridNTiles<decltype(gridDim)>(nFull, po.nTilesAssign);

	  std::cout << "gridDim: " << gridDim << " nTiles: " << nTiles << " nFull: " << nFull << std::endl;
	  const size_t nPixelsTile = helper::ii2n(tileDim);

	  hassertm3(data.size()==nFull*nPixelsTile, data.size(), nFull, nPixelsTile);
	  tileData.insert(tileData.end(), data.begin(), data.begin()+/*nTiles*/nFull*nPixelsTile);
	  nTiles=nFull;
	}
      else
	{
	  //
	  // synthetic data generation for testing
	  //



	  std::cout << __func__ << " generate synthetic data\n";

	  tileData.reserve(helper::ii2n(gridDim)*helper::ii2n(tileDim));

	  for(const auto gridId : helper::range2_n(gridDim))
	    for(const auto tileElemId : helper::range2_n(tileDim))
	      tileData.push_back(genSynthDataElem<D>(gridId, gridDim,
						     tileElemId, tileDim));

	  nTiles=helper::ii2n(gridDim);
	}
      return std::make_tuple(tileData, nTiles, tileDim, gridDim);
    }

    template<typename D, typename PO>
    auto initData(const PO& po)
    {
      return initData<D>(po,
			 po.configFNames.empty()
			 ? std::string("")
			 : po.configFNames.front());
    }

    template<typename I=uint32_t, typename DIM>
    auto read_qtLeafAssignment(const std::string loadAssignments, const DIM gridDim)
    {
      size_t nTiles=-1;
      std::vector<I> qtLeafAssignment;
      std::cout << "load assignment from file |" << loadAssignments << "|\n";
      //helper::readFile2(qtLeafAssignment, loadAssignments);
      //helper::bzip_decompress(qtLeafAssignment, loadAssignments);
      helper::readFileAuto(qtLeafAssignment, loadAssignments);
      std::cout << "loaded " << qtLeafAssignment.size() << " elements\n";
      if(qtLeafAssignment.empty())
	{
	  std::cout << "WARNING: ASSIGNMENT FILE "
		    << loadAssignments
		    << " NOT FOUND, CONTINUING WITH STANDARD INITIALIZATION\n\\n\n\n\n\n\n\n\n";
	}
      else if(nTiles!=static_cast<decltype(nTiles)>(-1)
	      && qtLeafAssignment.size() != nTiles)
	{
	  std::cerr << "ERROR: expected " <<  nTiles << " elements, read " << qtLeafAssignment.size() << std::endl;
	  throw std::runtime_error("error in loading file, wrong numbers of entries");
	}
      else if(helper::is_txt(loadAssignments))
	{
	  // FOR THE SAKE OF COMPARISON WITH OTHER APPROACHES,
	  // THE TXT FORMAT IS A LITTLE DIFFERENT
	  const auto grid2leaves=gridPos2QTLeaves(gridDim);	  
	  const auto data_txt=qtLeafAssignment;
	  for(const auto & i : helper::range_n(data_txt.size()))
	    qtLeafAssignment[grid2leaves[i]]=data_txt[i];
	}
      return qtLeafAssignment;
    }
  }
}

#endif //__SUPERTILES_PLACE_DATA__
