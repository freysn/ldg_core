#ifndef __SUPERTILES_PLACE_DRAW_TILE__
#define __SUPERTILES_PLACE_DRAW_TILE__

#include "helper_volData/vec.h"
#include "helper_CairoDraw.h"
#include "helper_imath.h"
#include "helper_idx.h"
#include "supertiles_configTypes.h"
#include "supertiles_QuadTree.h"
#include "helper_random.h"
#include "supertiles_place_cost.h"
#include "helper_color/cm_plasma.h"
#include "helper_color/cm_turbo.h"
#include "helper_color/cm_map.h"
#include "helper_blur.h"

#include "helper_pixelArtScaling.h"
#include "supertiles_place_DrawOpts.h"

#include "supertiles_place_useCases.h"



namespace supertiles
{
  namespace place
  {    
    
    template<typename AS, typename FTD, typename RTD, typename DIM, typename DIM_TILE>
    class DrawTile
    {

    public:

      using nodeId_t=uint32_t;

    private:
      using QT=supertiles_QuadTree<nodeId_t>;

      std::vector<uint8_t> _labels;

    public:

    DrawTile(const DIM gridDim_in,
	     const AS& qtLeafAssignment_in,
	     const FTD& feat_tileData_in,
	     const size_t nElemsFeatTile_in,
	     const RTD& rep_tileData_in,
	     const DIM_TILE rep_tileDim_in,
	     const repAggregationType_t repAggregationType_in,
	     const distFuncType_t distFuncType_in,
	     const std::vector<double>& disparities_in,
	     const DrawOpts& opts_in) :
      qt(gridDim_in),
      qtLeafAssignment(qtLeafAssignment_in),
      feat_tileData(feat_tileData_in),
      nElemsFeatTile(nElemsFeatTile_in),
      rep_tileData(rep_tileData_in),
      rep_tileDataDim(rep_tileDim_in),
      repAggregationType(repAggregationType_in),
      nodes2leaves(genMap_qt2leaves(qt)),
      distFuncType(distFuncType_in),
      disparities(disparities_in),
      opts(opts_in)
    {
      if(opts.mcmc_cm == 14831)
	{	  
	  const bool success=helper::bzip_decompress(_labels, opts.mcmc_labelsFName);
	  assert(success);
	}
    }

      auto tileOffsetsDims() const
      {
	const nodeId_t n=qt.nElems();
	std::vector<nodeId_t> offs(n+1);
	std::vector<DIM_TILE> dims(n);

	nodeId_t o=0;
	for(nodeId_t nodeId=0; nodeId<n; nodeId++)
	  {
	    dims[nodeId]=getTileDim(nodeId);
	    offs[nodeId]=o;
	    o += helper::ii2n(dims[nodeId]);
	  }
	offs[n]=o;

	return std::make_tuple(offs ,dims);
      }

      template<typename RO>
      static auto presentationVoidRepTileData(const size_t nElemsRepTile)
    {
      return std::vector<RO>(nElemsRepTile, RO());
    }


      DIM_TILE getTileDim(const nodeId_t nodeId) const
      {
	const auto level=qt.getLevel(nodeId);
	const auto scaleFac = helper::ipow2(level);

	DIM_TILE o;
	switch(repAggregationType)
	  {
	  case repAggregationType_chart:
	    {
	      // DIM_TILE dim(rep_tileDataDim.x, rep_tileDataDim.x);
	      // hassertm(rep_tileDataDim.y==1, rep_tileDataDim);

	      /*
	       used for stock1k
	      */
	      //DIM_TILE dim(32, 32);

	      //DIM_TILE dim(3, 3);
	      DIM_TILE dim(8, 8);

	      // const auto nLeaves=qt.nLeaves();

	      // if(nLeaves > 4*1024*1024)
	      // 	dim=DIM_TILE(3, 3);
	      o
		=
		dim*std::min(scaleFac, uint32_t(128));	      
	    }
	    break;
	  case repAggregationType_mcmc:
	    {
	      o=
		//max(rep_tileDataDim*scaleFac, decltype(rep_tileDataDim)(1,1));
		minv(maxv(decltype(rep_tileDataDim)(1,1)*scaleFac,
			  rep_tileDataDim)
		     , rep_tileDataDim*8);
	      break;
	    }
	  default:
	    o=rep_tileDataDim;
	  }
	return maxv(DIM_TILE(V2<double>(o)*opts.tileRepDimScale), DIM_TILE(8,8));
      }

      bool isLeafVoid(nodeId_t leafId) const
      {
	const auto tileId = qtLeafAssignment[leafId];
	hassertm3(tileId==voidTileIdx || tileId < rep_tileData.size(), leafId, tileId, rep_tileData.size());

	return (tileId == voidTileIdx);
      }

      std::vector<nodeId_t> getNonVoidLeaves(const std::vector<nodeId_t>& leaves) const
      {
	std::vector<nodeId_t> leavesSubset;
	leavesSubset.reserve(leaves.size());
	// leaf: leaves.size()==1
	for(const auto leafId : leaves)
	  if(!isLeafVoid(leafId))
	    leavesSubset.push_back(leafId);
	return leavesSubset;
      }

      std::vector<nodeId_t> getNonVoidLeaves(const nodeId_t nodeId) const
      {
	assert(nodeId<nodes2leaves.size());
	return getNonVoidLeaves(nodes2leaves[nodeId]);
      }

      template<typename IT>
      std::tuple<DIM_TILE, bool> operator()(IT to,
					    const nodeId_t nodeId,
					    const size_t shownNodesIdx=-1) const
      {

	
	
	using RO=
	  typename std::remove_const<typename std::remove_reference<decltype(to[0])>::type>::type;

	using ROE=
	  typename std::remove_const<typename std::remove_reference<decltype(to[0].x)>::type>::type;

	const auto rep_tileDimOut=getTileDim(nodeId);
	
	const auto & leaves = nodes2leaves[nodeId];


	//if(isLeaf && isLeafVoid(nodeId))
	if(getNonVoidLeaves(nodes2leaves[nodeId]).empty())
	  {
	    const auto nElemsRepTileOut=helper::ii2n(rep_tileDimOut);
	    const auto voidRepTileData =
	      presentationVoidRepTileData<RO>(nElemsRepTileOut);
	    auto from_it = voidRepTileData.begin();
	    std::copy(from_it,
		      from_it+nElemsRepTileOut,
		      to
		      );
	    return std::make_tuple(rep_tileDimOut, false);

	  }

	auto tileIdFront = [&]()
	{
	  for(const auto & leafId : leaves)
	    {
	      if(isLeafVoid(leafId))
		continue;
	      const auto tileId = qtLeafAssignment[leafId];
	      assert(tileId != voidTileIdx);
	      return tileId;
	    }
	  return voidTileIdx;
	};

	auto caltechCopy = [&](nodeId_t tileId, double disparity)
	{
	  if(tileId != voidTileIdx)
	    {
	      
	      const auto nElemsRepTile=helper::ii2n(rep_tileDataDim);

	      using SRC=
		typename std::remove_const<typename std::remove_reference<decltype(rep_tileData[0])>::type>::type;

	      const SRC* src=&rep_tileData[tileId*nElemsRepTile];

	      std::vector<SRC> tmp;
	      tmp.reserve(helper::ii2n(rep_tileDimOut));
	      
	      if(rep_tileDataDim!=rep_tileDimOut)
		{		  
		  DIM_TILE p;
		  for(p.y=0; p.y<rep_tileDimOut.y; p.y++)
		    for(p.x=0; p.x<rep_tileDimOut.x; p.x++)
		    {
		      V2<int> a=
			(V2<double>(p.x+.5, p.y+.5)
			 /V2<double>(rep_tileDimOut))
			*V2<double>(rep_tileDataDim);

		      hassertm(a.x>=0 && a.x < rep_tileDataDim.x
				&&
				a.y>=0 && a.y < rep_tileDataDim.y
				,a);
		      tmp.push_back(src[a.x+rep_tileDataDim.x*a.y]);		      
		    }
		  src=&tmp[0];
		}
	      
	      struct Op
	      {
		void operator()(V4<double>* dest,
				const SRC* src,
				const size_t nElemsRepTile) const
		{
		  const auto buf =
		    //&rep_tileData[tileId*nElemsRepTile];
		    caltech2col
		    (src,
		     nElemsRepTile);

		  std::copy(buf.begin(), buf.end(), dest);
		}

		void operator()(V4<uint8_t>* dest,
				const V4<uint8_t>* src,
				const size_t nElemsRepTile) const
		{
		  std::copy(&src[0], &src[nElemsRepTile], dest);
		}

		void operator()(V4<uint8_t>* /*dest*/,
				const double* /*src*/,
				const size_t /*nElemsRepTile*/) const
		{
		  // wrong code path
		  assert(false);
		}

		void operator()(double* /*dest*/,
				const double* /*src*/,
				const size_t /*nElemsRepTile*/) const
		{
		  // wrong code path
		  assert(false);
		}
	      };

	      Op()(&to[0], src,
		   helper::ii2n(rep_tileDimOut));

	      // std::cout << "BLUR nodeId " << nodeId
	      // 		<< " disparity " << disparity
	      // 		<< std::endl;

	      if(opts.blur && disparity > 0.)
		{
		  hassertm(disparity <= 1., disparity);

		  const auto radius=
		    std::min(30.,
		    std::max(rep_tileDimOut.x, rep_tileDimOut.y)
			     *0.05
			     *disparity);

		  std::cout << "BLUR nodeId " << nodeId
			    << " disparity " << disparity
			    << " blurRadius " << radius
			    << std::endl;

		  helper::BlurFunctor blur(rep_tileDimOut.x,
					   rep_tileDimOut.y,
					   radius);


		  const std::vector<RO>buf_orig(to, to+nElemsRepTile);
#ifndef NO_OMP
#pragma omp parallel for
#endif
		  for(size_t pxIdx=0; pxIdx<nElemsRepTile; pxIdx++)
		    to[pxIdx] = blur(pxIdx, buf_orig);
		}


	      // for(const auto & i : helper::range_n(nElemsRepTile))
	      //   to[i] = buf[i];
	    }
	};


	auto getCentralTileId = [&]()
	{
	  auto leavesSubset=getNonVoidLeaves(nodes2leaves[nodeId]);

	  const auto nonVoidCnt=leavesSubset.size();
	  if(nonVoidCnt==1 || nonVoidCnt==2)
	    {
	      const auto tileId = tileIdFront();
	      assert(tileId==qtLeafAssignment[leavesSubset.front()]);
	      return tileId;
	    }
	  else if(nonVoidCnt>2)
	    {
	      if(true)
		{

		  const size_t maxNLeaves=1000;
		  if(leavesSubset.size()>maxNLeaves)
		    {
		      //std::cout << "resize " << leavesSubset.size() << " to " << maxNLeaves << std::endl;
		      std::mt19937 rng(1337);
		      helper::mrandom_shuffle(leavesSubset.begin(), leavesSubset.end(), rng);
		      leavesSubset.resize(maxNLeaves);
		    }
		}

	      const auto up = helper::uniquePairs(leavesSubset.size());
	      std::vector<double> dists(up.size(), 0.);

	      // if(up.size()>10000)
	      // 	std::cout << "#unique pairs: " << up.size() << "(#leaves: " << leavesSubset.size() << ")" << std::endl;

#pragma omp parallel for
	      for(size_t pairIdx=0; pairIdx<up.size(); pairIdx++)
		{
		  const auto lp=std::make_pair(leavesSubset[up[pairIdx].first],
					       leavesSubset[up[pairIdx].second]);

		  assert(!(isLeafVoid(lp.first) || isLeafVoid(lp.second)));

		  const auto lpt=
		    std::make_pair(qtLeafAssignment[lp.first],
				   qtLeafAssignment[lp.second]);

		  //auto from=rep_tileData.begin()+l*nElemsRepTile;


		  auto from=
		    std::make_pair
		    // (feat_tileData.begin()+lpt.first*nElemsFeatTile,
		    //  feat_tileData.begin()+lpt.second*nElemsFeatTile);
		    (get_feat_tileData(lpt.first), get_feat_tileData(lpt.second));

		  // auto from=
		  //   std::make_pair
		  //   (feat_supertileData.begin()+lp.first*nElemsFeatTile,
		  //    feat_supertileData.begin()+lp.second*nElemsFeatTile);

		  //DistOp<distFuncType,double> distOp;

		  auto applyDistOp = [&](auto distOp)
		  {
		    for(const auto & i : helper::range_n(nElemsFeatTile))
		      distOp(from.first[i],from.second[i]);
		    dists[pairIdx]=distOp.get();
		  };

		  switch(distFuncType)
		    {
		    case distFuncType_cosine_normalized:
		      applyDistOp(DistOp<distFuncType_cosine_normalized,double>());
		      break;
		    case distFuncType_norm2:
		      applyDistOp(DistOp<distFuncType_norm2,double>());
		      break;
		    default:
		      assert(false);
		      break;
		    };
		}

	      //
	      // pick leaf with least distance to others
	      //
	      std::pair<double, size_t>
		best(std::numeric_limits<double>::max(),
		     -1);


	      for(const auto l0 : helper::range_n(leavesSubset.size()))
		{
		  if(isLeafVoid(leavesSubset[l0]))
		    continue;
		  double d=0.;
		  for(const auto l1 : helper::range_n(leavesSubset.size()))
		    if(l0!=l1)
		      d+=dists[helper::uniquePairIdx(std::min(l0,l1),
						     std::max(l0,l1),
						     leavesSubset.size())];
		  if(d<best.first)
		    best=std::make_pair(d, l0);
		}

	      hassertm2(best.second!=static_cast<size_t>(-1), best.first, best.second);
	      const auto leafId = leavesSubset[best.second];
	      hassertm2(!isLeafVoid(leafId), best.first, best.second);
	      const auto tileId = qtLeafAssignment[leafId];
	      hassertm2(tileId!=voidTileIdx, best.first, best.second);
	      return tileId;
	    }
	  return voidTileIdx;
	};


	if(repAggregationType==repAggregationType_front)
	  {
	    // just use value of first leaf
	    //throw "repAggregationType_front not supported at the moment";

	    assert(nodeId<disparities.size());

	    const auto disparity=disparities[nodeId];
	    const auto tileId = tileIdFront();
	    caltechCopy(tileId, disparity);

	  }
	else
	  if(repAggregationType==repAggregationType_average_colRGB ||
	     repAggregationType==repAggregationType_mcmc)
	    {
	      const bool mcmc=(repAggregationType==repAggregationType_mcmc);
	      const bool colRGB=(repAggregationType==repAggregationType_average_colRGB);

	      const auto nElemsRepTile=helper::ii2n(rep_tileDataDim);
	      RTD buf(nElemsRepTile, typename RTD::value_type());

	      //static_assert(std::is_same<typename RTD::value_type, double>::value, "should be double");

	      // counts non void tiles involved
	      size_t cnt=0;
	      for(const auto & li : leaves)
		{
		  if(isLeafVoid(li))
		    continue;
		  const auto l=qtLeafAssignment[li];
		  cnt++;
		  auto from=rep_tileData.begin()+l*nElemsRepTile;
		  for(const auto & i : helper::range_n(nElemsRepTile))
		    {
		      auto e = from[i];
		      if(mcmc)
			e = mcmc_binarize(e);
		      buf[i]+=e;
		    }
		}

	      if(cnt > 0)
		{
		  for(const auto & i : helper::range_n(nElemsRepTile))
		    buf[i]/= cnt;

		  uint32_t ifac =
		    helper::nextPow2(
				     std::max(
					      rep_tileDimOut.x/rep_tileDataDim.x,
					      rep_tileDimOut.y/rep_tileDataDim.y));

		  const auto ifac_pow2 = helper::ilog2(ifac);
		  //hassertm(helper::isPow2(ifac), ifac);

#ifdef USE_STRAX_MCMC
		  const bool doVR25D=false;//(cnt>=4096)/*cnt>128*/;

		  if(doVR25D)
		    {
		      //
		      // INCORPORATE VOLUME RENDERING STUFF FROM STRAX HERE
		      //
		      strax_mcmc_opts<> opts;
		      opts.imgDim=rep_tileDimOut;
		      //opts.epxNIters=std::min(ifac_pow2, 1u);
		      opts.epxNIters=0;
		      //opts.epxNIters=1;
		      opts.tstepModifier=0.55;
		      opts.rko.nRaysPerPixel=32;
		      opts.rko.nRaysAmbientOcclusion=8;
		      //opts.rko.maxLenFacAmbientOcclusion=0.02;
		      opts.colMap.b=V4<double>(230/255., 0., 126/255., 1.);
		      opts.camOrigin=V3<double>(0.,-0.15,1.)*2.7;
		      std::random_device r;
		      std::mt19937 rng(r());
		      const auto rslt=strax_mcmc<false>(buf, rep_tileDataDim,
							opts, rng);

		      const auto nElemsRepTileOut=helper::ii2n(rep_tileDimOut);
		      //std::copy(&rslt[0], &rslt[nElemsRepTileOut], to);
		      
		      
		      for(const auto & i : helper::range_n(nElemsRepTileOut))
			to[i]=255.*rslt[i];

		      //std::cout << "rep_tileDimOut " << rep_tileDimOut << std::endl;
		    }
		  else
#endif // USE_STRAX_MCMC
		    {
		      if(rep_tileDataDim!=rep_tileDimOut)
			{
			  auto dim=rep_tileDataDim;
			  if(ifac > 0)
			    {
			      
			      //
			      // UPSAMPLE
			      //
			      // hassertm3(ifac*rep_tileDataDim
			      // 	    ==rep_tileDimOut,
			      // 	    ifac, rep_tileDataDim,
			      // 	    rep_tileDimOut);

			      for(size_t e=0; e<ifac_pow2; e++)
				std::tie(buf, dim)
				  =helper::epx(buf.begin(),
					       dim);
			      if(helper::ii2n(dim)>=256*256)
				std::cout << "WARNING: LARGE DIM" << dim << " " << ifac_pow2 << " " << rep_tileDataDim << " " << rep_tileDimOut<< std::endl;
			    }



			  // hassertm2(dim.x>=rep_tileDimOut.x &&
			  // 		dim.y>=rep_tileDimOut.y,
			  // 		dim, rep_tileDimOut);

			  {
			    helper::BlurFunctor bf
			      (dim.x, dim.y, std::max(ifac, static_cast<uint32_t>(1)));

			    RTD buf2=buf;

			    // for(const auto& i : helper::range_n(buf.size()))
			    //   buf2[i]=bf(i, buf);

			    buf.resize(helper::ii2n(rep_tileDimOut));
			    for(const auto & id : helper::range2_n(rep_tileDimOut))
			      {
				using D2=V2<double>;
				const auto np
				  =
				  D2(id)
				  /
				  D2(rep_tileDimOut);

				auto id_buf2 =
				  DIM_TILE(np*D2(dim));

				buf[helper::ii2i(id, rep_tileDimOut)]
				  =
				  bf(helper::ii2i(id_buf2, dim), buf2);
			      }

			  }
			}		    
		  
		      const auto nElemsRepTileOut=helper::ii2n(rep_tileDimOut);

		      hassertm2(buf.size()==nElemsRepTileOut,
				buf.size(), rep_tileDimOut);

		      struct CopyOp
		      {
			void operator()(V4<uint8_t>* dest,
					const V4<double>* src,
					const size_t nElemsRepTileOut) const
			{
			  for(const auto & i : helper::range_n(nElemsRepTileOut))
			    dest[i]=255.*src[i];
			}

			void operator()(V4<double>* dest,
					const V4<double>* src,
					const size_t nElemsRepTileOut) const
			{
			  std::copy(&src[0], &src[nElemsRepTileOut], dest);
			}

			void operator()(double* /*dest*/,
					const V4<double>* /*src*/,
					const size_t /*nElemsRepTileOut*/) const
			{
			  assert(false);
			}
		      };

		  if(mcmc)
		    {
		      if(opts.mcmc_cm != 14831)
			{
			  const auto rslt = mcmc2col(&buf[0], nElemsRepTileOut, opts.mcmc_cm);
			  
			  CopyOp()(&to[0], &rslt[0],
				   nElemsRepTileOut);
			}
		      else
			{
			  //
			  // SPECIAL MODE FOR COLORING LABELLED DATA
			  //

			  // from gen_glyph.py
			  const std::vector<V3<double>> labelColMap={{0.9019607843137255,0.09803921568627451,0.29411764705882354}, {0.23529411764705882,0.7058823529411765,0.29411764705882354}, {1.0,0.8823529411764706,0.09803921568627451}, {0.2627450980392157,0.38823529411764707,0.8470588235294118}, {0.9607843137254902,0.5098039215686274,0.19215686274509805}, {0.5686274509803921,0.11764705882352941,0.7058823529411765}, {0.25882352941176473,0.8313725490196079,0.9568627450980393}, {0.9411764705882353,0.19607843137254902,0.9019607843137255}, {0.7490196078431373,0.9372549019607843,0.27058823529411763}, {0.9803921568627451,0.7450980392156863,0.8313725490196079}, {0.27450980392156865,0.6,0.5647058823529412}, {0.8627450980392157,0.7450980392156863,1.0}, {0.6039215686274509,0.38823529411764707,0.1411764705882353}, {1.0,0.9803921568627451,0.7843137254901961}, {0.5019607843137255,0.0,0.0}, {0.6666666666666666,1.0,0.7647058823529411}, {0.5019607843137255,0.5019607843137255,0.0}, {1.0,0.8470588235294118,0.6941176470588235}, {0.0,0.0,0.4588235294117647}, {0.6627450980392157,0.6627450980392157,0.6627450980392157}};;
			  
			  assert(!_labels.empty());
			  std::map<uint8_t, size_t> labelCounts;
			  			  
			  
			  for(const auto & li : leaves)
			    {
			      if(isLeafVoid(li))
				continue;

			      const auto l=qtLeafAssignment[li];
			      
			      hassertm2(l<_labels.size(), l, _labels.size());
			      const auto o=_labels[l];
			      
			      auto it=labelCounts.find(o);
			      if(it != labelCounts.end())
				it->second++;
			      else
				labelCounts.emplace(o,1);
			    }

			  assert(!labelCounts.empty());
			  
			  size_t totalCnt=0;
			  {
			    size_t cnt=0;			    
			    for(auto & e : labelCounts)
			      {
				e.second+=cnt;
				cnt=e.second;
			      }
			    totalCnt=cnt;
			  }
			  			  
			  assert(totalCnt>0);
			  
			  V2<uint32_t> dim;
			  dim.x=std::sqrt(nElemsRepTileOut)+.5;
			  dim.y=dim.x;

			  assert(dim.x*dim.y==nElemsRepTileOut);

			  std::vector<V4<double>> rslt(dim.x*dim.y);


			  
			  auto it=labelCounts.begin();			  

			  const auto nLabeledMembersPerRow=static_cast<double>(totalCnt)/static_cast<double>(dim.y);
			  
			  for(const auto & y : helper::range_n(dim.y))
			    {
			      
			      while(it != labelCounts.end() && it->second < y*nLabeledMembersPerRow)
				it++;				

			      assert(it != labelCounts.end());
			      
			      const auto labelId=it->first;
			      
			      const V3<double> labelCol=labelColMap[labelId];
			      const V3<double> baseCol=labelColMap.back();

			      
			      
			      for(const auto & x : helper::range_n(dim.x))
				{
				  const size_t idx=y*dim.x+x;

				  struct GetCol
				  {

				    GetCol(V3<double> labelColIn, V3<double> baseColIn)
				      :
				      labelCol(labelColIn),
				      baseCol(baseColIn)
				    {};
				    
				    auto operator()(double b)
				    {
				      return b*labelCol + (1-b)*baseCol;
				    };

				    auto operator()(V4<uint8_t> b)
				    {
				      return (b.x*labelCol + (255-b.x)*baseCol)/255.;
				    }

				    auto operator()(V3<float> b)
				    {
				      return this->operator()(b.x);
				    }

				    const V3<double> labelCol;
				    const V3<double> baseCol;				    
				  };
				  
				  const V3<double> col= GetCol(labelCol, baseCol)(buf[idx]);
				    
				  rslt[idx]=V4<double>(col.x, col.y, col.z, 1.);
				}
			    }
			  
			  
			  CopyOp()(&to[0], &rslt[0],
				   nElemsRepTileOut);
			  
			}
		    }
		  else if(colRGB)
		    {
		      // this should be the rgb case
		      // currently not supported, different number of input and output elements
		      // (3 to 1)

		      //assert(nElemsRepTileOut==1);

		      const auto rslt = colRGB2col(&buf[0], nElemsRepTileOut);
		      CopyOp()(&to[0], &rslt[0],
			       nElemsRepTileOut);

		      // for(const auto & i : helper::range_n(nElemsRepTileOut))
		      // 	{
		      // 	  to[i]=rslt[i];
		      // 	}
		    }
		  else
		    assert(false);
		}
		}
	      else
		{
		  // copy empty tile
		  assert(isLeafVoid(leaves.front()));
		  //copyTile(leaves.front());
		}
	    }
	  else if (repAggregationType==repAggregationType_chart)
	    {
	      using F2=V2<int>;
	      using CairoDraw_t = helper::CairoDraw<F2>;

	      F2 chartBoxDim(rep_tileDimOut);

	      
	      std::vector<nodeId_t> tileIds;
	      
	      {
		const auto leaves=getNonVoidLeaves(nodes2leaves[nodeId]);
		tileIds.reserve(leaves.size());
		for(const auto e : leaves)
		  tileIds.push_back(qtLeafAssignment[e]);
	      }

	      CairoDraw_t cd(helper::cairoBackend_array, chartBoxDim);
	      auto cr=cd.get();
	      drawSignalAll(cr, tileIds, feat_tileData, nElemsFeatTile,
			      chartBoxDim, cm_turbo, opts);
		

		if(false)
		  {
		    //const auto colMap=cm_plasma;
		    const auto centralTileId=getCentralTileId();

		    if(centralTileId != voidTileIdx)
		      drawSignal(cr, feat_tileData.begin()+centralTileId*nElemsFeatTile, nElemsFeatTile, chartBoxDim, cm_turbo, opts);
		  }
	      const auto [img,dim] = cd.getDataRGBA<ROE>();
	      hassertm2(dim==chartBoxDim, dim, chartBoxDim);

	      std::copy(img.begin(), img.end(), to);

	    }
	  else if(repAggregationType==repAggregationType_central)
	    {
	      assert(nodeId<disparities.size());

	      const auto disparity=disparities[nodeId];
	      caltechCopy(getCentralTileId(), disparity);
	    }
	  else
	    assert(false);
	
	    
	if(shownNodesIdx<opts.shownNodes.size())
	  {
	    const auto markMode= opts.shownNodes[shownNodesIdx][nodeId];
	    if((markMode==2/*|| opts.shownNodes[nodeId]==3*/) &&
	       opts.shownNodes_inactiveOpMod<1.)
	      {
		for(auto it=to; it != to+helper::ii2n(rep_tileDimOut); it++)
		  it->w*=opts.shownNodes_inactiveOpMod;
	      }
	    else if(markMode==2/*|| opts.shownNodes[nodeId]==3*/ && false)
	      {
		for(auto it=to; it != to+helper::ii2n(rep_tileDimOut); it++)
		  {
		    typename
		      std::remove_reference<decltype(*it)>::type tmp(0.6, 0.3,1., 0.4);
		    std::swap(tmp, *it);
		    over(*it,tmp);
		  }
	      }
	    else if(markMode==3 && false)
	      {
		for(auto it=to; it != to+helper::ii2n(rep_tileDimOut); it++)
		  {
		    typename
		      std::remove_reference<decltype(*it)>::type tmp(0.3,1.,0.6, 0.4);
		    std::swap(tmp, *it);
		    over(*it, tmp);
		  }
		/*
		auto cr=cd.get();

		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		//cairo_set_line_width(cr, tileImg.lineWidth);
		cairo_set_line_width(cr, 10);

		V2<double> imgPos, imgDim;
		std::tie(imgPos,
			 imgDim,
			 std::ignore,
			 std::ignore)
		  = node2tileImg_(nodeId, false);


		cairo_rectangle(cr, p.x, p.y+dim.y, f*dim.x, dim.y);

		cairo_stroke(cr);
		*/
		//   }	      
	      }
	  }
      
	return std::make_tuple(rep_tileDimOut, true);
      }

      auto getRepAggregationType() const
      {
	return repAggregationType;
      }

      auto get_feat_tileData(nodeId_t tileId) const
      {
	hassertm2(tileId*nElemsFeatTile < feat_tileData.size(), tileId,  feat_tileData.size()/nElemsFeatTile);
	return feat_tileData.begin()+tileId*nElemsFeatTile;
      }

      auto get_feat_tileData_fromLeaf(nodeId_t leafId) const
      {
	assert(!isLeafVoid(leafId));
	const auto tileId=qtLeafAssignment[leafId];
	return get_feat_tileData(tileId);
      }

      auto get_nElemsFeatTile() const
      {
	return nElemsFeatTile;
      }

      const std::vector<nodeId_t>& node2leaf(nodeId_t nodeId) const
      {
	assert(nodeId<nodes2leaves.size());
	return nodes2leaves[nodeId];
      }
      
    private:
      const QT qt;

      const AS& qtLeafAssignment;
      const FTD& feat_tileData;
      const size_t nElemsFeatTile;
      const RTD& rep_tileData;
      const DIM_TILE rep_tileDataDim;
      const repAggregationType_t repAggregationType;

      std::vector<std::vector<nodeId_t>> nodes2leaves;

      const distFuncType_t distFuncType;
      const std::vector<double>& disparities;

      const DrawOpts& opts;
    };
  }
}

#endif // __SUPERTILES_PLACE_DRAW_TILE__
