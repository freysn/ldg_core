#ifndef __SUPERTILES_PLACE_PLAN__
#define __SUPERTILES_PLACE_PLAN__

#include "helper_random.h"

namespace supertiles
{
  namespace place
  {


    
    class Plan
    {
      enum class planMode_t
	{
	  simpleAll, inBoxes, acrossNodeHierarchies
	};

      static const planMode_t planMode=planMode_t::simpleAll;
      
      using IDX=uint32_t;
      using DIM=V2<IDX>;
    public:
      Plan(uint32_t level_in,
	   bool aggregateExchangeMode_in,
	   const std::vector<uint32_t>& planIdcs_in,
	   const std::vector<IDX>& grid2qt
	   ) :
	level(level_in),
	aggregateExchangeMode(aggregateExchangeMode_in),
	planIdcs(planIdcs_in)
      {	
	if(planMode==planMode_t::inBoxes)
	  {

	    //"TODO: NEED TO HANDLE NON-POW4 GRIDS"
	    assert(helper::ilog4_ceil(planIdcs.size())
		   ==helper::ilog4_floor(planIdcs.size()));
	    
	    const size_t nLevels=helper::ilog4_ceil(planIdcs_in.size());

	    auto prob =[this](size_t level)
	    {
	      size_t start=0;

	      if(!inBoxes_planProbList.empty())
		start=std::get<2>(inBoxes_planProbList.back());
	      return start+(1<<level);
	    };
	    //const auto nLeaves = grid2qt.size();
	    
	    const size_t nLeaves=helper::ipow4(nLevels);
	
	    //assert(nLevels==helper::ilog4_floor(nLeaves));

	    const DIM gridDim(helper::ipow2(nLevels),
			      helper::ipow2(nLevels));

	
	    DIM grid2qtDim;
	    {
	      const size_t l=helper::ilog4_ceil(grid2qt.size());
	      grid2qtDim.x=helper::ipow2(l);
	      grid2qtDim.y=helper::ipow2(l);
	    }

	    //assert(helper::ii2n(gridDim)==grid2qt.size());


	    std::cout << "???" << nLeaves << " " << gridDim << " " << planIdcs_in.size() << " " << grid2qtDim << std::endl;


	
	
	
	    for(const auto planLevel : helper::range_be(1,nLevels))
	      for(const bool boxShiftY : {false, true})
		for(const bool boxShiftX : {false, true})
		  {
		    const V2<bool> boxShift(boxShiftX, boxShiftY);

		    std::string idStr;
		    {
		      std::stringstream ss;
		      ss<<"(relative) level " << planLevel
			<< ", boxShift " << boxShift;
		      idStr=ss.str();
		    }		
		    std::cout << "\n\n\n\n\n\n\n\n\n\n############### generate plan at " << idStr << std::endl;
		    DIM boxDim;
		    {
		  
		      boxDim.x=helper::ipow2(planLevel);
		      boxDim.y=boxDim.x;
		    }

		

		    std::vector<IDX> plan;
		    std::vector<IDX> offsets;		

		    auto finishPlanGroup= [&plan, &offsets]()
		    {		  
		      offsets.push_back(plan.size());
		      //std::cout << "-----------finish plan group at index " << offsets.back() << idStr << std::endl;
		    };

		    auto addNodePosIfValid= [&plan, &grid2qt, &gridDim, &grid2qtDim](auto p)
		    {
		      //std::cout << "addPos " << p << std::endl;
		      hassertm4(p.x<gridDim.x && p.y < gridDim.y, p.x, p.y, gridDim.x, gridDim.y);

		      const auto gridIdx=helper::ii2i(p, /*gridDim*/grid2qtDim);
		      const auto nodeId=grid2qt[gridIdx];
		      
		      
		      //std::cout << "add nodeId " << nodeId << " at index " << plan.size() << std::endl;
		      plan.push_back(nodeId);
		    };

		    std::cout << "========== ADD PLAN BODY ==========\n";
		    {
		      const V2<IDX> shiftGrid((boxShift.x ? boxDim.x/2 : 0),
					      (boxShift.y ? boxDim.y/2 : 0));
		      DIM off;	
		      for(off.y=shiftGrid.y; off.y<gridDim.y-shiftGrid.y; off.y+=boxDim.y)
			for(off.x=shiftGrid.x; off.x<gridDim.x-shiftGrid.x; off.x+=boxDim.x)
			  {
			    for(const auto & e : helper::range2_n<IDX>(boxDim))
			      {
				const DIM p(off+e);
				//std::cout << "gridPos: " << p << " " << off << " " << e << std::endl;
				addNodePosIfValid(p);		
			      }
			    finishPlanGroup();
			  }       
		    }

		    auto addBordersHorizontal= [addNodePosIfValid, finishPlanGroup](auto gridDim, auto boxDim, auto boxShift, bool swapXY)
		    {
		      if(swapXY)
			{
			  std::swap(gridDim.x, gridDim.y);
			  std::swap(boxDim.x, boxDim.y);
			  std::swap(boxShift.x, boxShift.y);
			}
		      if(boxShift.x)
			{
			  const auto off_y=boxShift.y*(boxDim.y/2);
			  for(auto off : {DIM(0,off_y), DIM(gridDim.x-boxDim.x/2, off_y)})
			    for(; off.y<gridDim.y-off_y; off.y+=boxDim.y)
			      {
				for(const auto & e : helper::range2_n<IDX>(DIM(boxDim.x/2, boxDim.y)))
				  {
				    auto p = off+e;
				    if(swapXY)
				      std::swap(p.x, p.y);
				    addNodePosIfValid(p);
				  }
				finishPlanGroup();
			      }
			}
		    };

		
		    for(const bool swapXY : {false, true})
		      {
			std::cout << "========== ADD BORDERS HORIZONTAL ========== " << swapXY << " " << idStr << "\n";
			addBordersHorizontal(gridDim, boxDim, boxShift, swapXY);
		      }


		    std::cout << "========== ADD BOX CORNERS ========== " << idStr << "\n";
		    // add box corners
		    if(boxShift.x&&boxShift.y)
		      {
			for(const auto off : {DIM(0,0),
					      DIM(gridDim.x-boxDim.x/2, 0),
					      DIM(0, gridDim.y-boxDim.y/2),
					      DIM(gridDim.x-boxDim.x/2, gridDim.y-boxDim.y/2)})
			  {
			    for(const auto & e : helper::range2_n<IDX>(DIM(boxDim.x/2, boxDim.y/2)))
			      addNodePosIfValid(off+e);
			    finishPlanGroup();
			  }
		      }
		
		
		    hassertm2(plan.size()==planIdcs.size(), plan.size(), planIdcs.size());

#ifndef NDEBUG
		    {
		      auto cpyA=planIdcs;
		      auto cpyB=plan;
		      std::sort(cpyA.begin(), cpyA.end());
		      assert(cpyA==helper::range_n<std::vector<IDX>>(nLeaves));
		      std::sort(cpyB.begin(), cpyB.end());
		      // for(const auto i : helper::range_n(cpyB.size()))
		      //   std::cout << i << ": " << cpyB[i] << std::endl;
		      assert(cpyB==helper::range_n<std::vector<IDX>>(nLeaves));
		      assert(cpyA==cpyB);
		    }
#endif		
		    inBoxes_planProbList.emplace_back(plan, offsets, prob(planLevel));
		    // for(const auto e : offsets)
		    //   std::cout << "OFFSET " << e << std::endl;
		  }
	    inBoxes_planProbList.emplace_back(planIdcs, std::vector<IDX>(1, planIdcs.size()), prob(nLevels));
	  }
	
      }
  
      template<typename RNG>
      auto operator()(RNG& rng)
      {
	size_t leafLevel=0;	

	
	
	if(aggregateExchangeMode)
	  {	    
	    leafLevel=level;
	  }


	switch(planMode)
	  {
	  case planMode_t::acrossNodeHierarchies:
	    {
	const size_t nLeaves=planIdcs.size();
	const size_t nLevels=helper::ilog4_ceil(nLeaves);

	std::cout << "nLeaves: " << nLeaves << "; nLevels: " << nLevels << std::endl;
	std::vector<std::tuple<size_t, size_t>> stridesProbs;
	size_t probsSum=0;
	for(size_t level=0; level<nLevels; level++)
	  {
	    std::cout << "stride: " << helper::ipow4(level) << std::endl;

	    const size_t prob=4;//(1<<level);
	    probsSum+=prob;
	    stridesProbs.emplace_back(helper::ipow4(level), probsSum);	    
	  }

	size_t stride=0;
	{
	  std::uniform_int_distribution<IDX> dis(0, probsSum);
	  const auto sel=dis(rng);
	  size_t idx=0;
	  while(sel>std::get<1>(stridesProbs[idx]))
	    idx++;
	  assert(idx < stridesProbs.size());
	  stride=std::get<0>(stridesProbs[idx]);
	}


	std::cout << "SELECTED STRIDE: " << stride << std::endl;
	auto tmp = helper::range_n<std::vector<IDX>>(helper::ipow4(nLevels));
	
	for(size_t i=0; i<tmp.size(); i+=stride)
	  {
	    helper::mrandom_shuffle(tmp.begin()+i, tmp.begin()+i+stride, rng);
	  }
	

	std::fill(planIdcs.begin(), planIdcs.end(), -1);
	const size_t partitionSize=(4*stride);
	const size_t nPartitions=nLeaves/partitionSize;
	
	assert(nPartitions>0);
	for(size_t i=0; i<nPartitions; i++)
	  {
	    const auto baseOff=i*partitionSize;
	    for(size_t j=0; j<stride; j++)
	      for(size_t k=0; k<4; k++)
		{
		  const auto offTo=4*j+k;
		  const auto offFrom=j+k*stride;
		  
		  const auto idxTo=baseOff+offTo;
		  const auto idxFrom=baseOff+offFrom;
		  //std::cout << idxTo << " <= " << idxFrom << " " << tmp[idxFrom] << "(stride: " << stride << "; nLeaves: " << nLeaves << ")"<<std::endl;
		  hassertm2(offTo<partitionSize, offTo, partitionSize);
		  hassertm2(offFrom<partitionSize, offFrom, partitionSize);
		  planIdcs[idxTo]=tmp[idxFrom];
		}
	  }
#ifndef NDEBUG
	{
	  auto tmp = helper::range_n<std::vector<IDX>>(helper::ipow4(nLevels));
	  auto tmp2=planIdcs;
	  std::sort(tmp2.begin(), tmp2.end());

	  // for(const auto i : helper::range_n(tmp2.size()))
	  //   std::cout << "i " << i << ": " << tmp2[i] << std::endl;
	  assert(tmp==tmp2);
	}
#endif
	// for(const auto i : helper::range_n(planIdcs.size()))
	//   std::cout << "i " << i << ": " << planIdcs[i] << std::endl;
	return std::make_tuple(planIdcs.cbegin(), planIdcs.cend(),
			       aggregateExchangeMode, leafLevel);
	    }
	  case planMode_t::inBoxes:
	    {

	//"USE SPATIAL SIMILARITY AS CRITERION, ADDITIONALLY REQUIRE GRID POS MAP"


	assert(!inBoxes_planProbList.empty());
	std::uniform_int_distribution<IDX> dis(0,std::get<2>(inBoxes_planProbList.back()));

	size_t planListIdx=0;
	{	  
	  const auto sel=dis(rng);
	  while(sel>std::get<2>(inBoxes_planProbList[planListIdx]))
	    planListIdx++;
	  assert(planListIdx < inBoxes_planProbList.size());	  
	}
	std::cout << "planListIdx " << planListIdx << " (of " << inBoxes_planProbList.size() << ")" << std::endl;
	
	auto & plan = std::get<0>(inBoxes_planProbList[planListIdx]);
	auto & offsets = std::get<1>(inBoxes_planProbList[planListIdx]);
	

	for(size_t i=0; i<offsets.size(); i++)
	  {
	    IDX b=0;
	    if(i>0)
	      b=offsets[i-1];

	    const auto e = offsets[i];
	    assert(b<=e);
	    assert(e<=plan.size());

	    helper::mrandom_shuffle(plan.begin()+b, plan.begin()+e, rng);
	  }

	return std::make_tuple(plan.cbegin(), plan.cend(),
			       aggregateExchangeMode, leafLevel);
	    }

	// const auto stride=helper::ipow4(leafLevel);
// 	if(stride>1)
// 	  {
// #pragma omp parallel for
// 	    for(size_t i=0; i<costsLeaves.size();i+=stride)
// 	      {
// 		auto & o = costsLeaves[i];
// 		for(size_t j=1; j<stride; j++)
// 		  o+=costsLeaves[i+j];
// 	      }
// 	  }

// 	const size_t n=0.15*planIdcs.size();
// 	assert(n<planIdcs.size()/2);

// 	// create two partitions of elements with large and low costs
// 	std::partial_sort(planIdcs.begin(), planIdcs.begin()+n, planIdcs.end(),
// 			  [&](auto i0, auto i1){return costsLeaves[i0*stride] > costsLeaves[i1*stride];});


// 	// shuffle elements with lower cost
// 	helper::mrandom_shuffle(planIdcs.begin()+n, planIdcs.end(), rng);

// 	// shuffle full plan
// 	helper::mrandom_shuffle(planIdcs.begin(), planIdcs.begin()+2*n, rng);
	  
// 	return std::make_tuple(planIdcs.cbegin(), planIdcs.cbegin()+2*n, aggregateExchangeMode, leafLevel);

	  case planMode_t::simpleAll:
	    helper::mrandom_shuffle(planIdcs.begin(), planIdcs.end(), rng);
	    return std::make_tuple(planIdcs.cbegin(), planIdcs.cend(),
				   aggregateExchangeMode, leafLevel);
	  default:
	    assert(false);
	    return std::make_tuple(planIdcs.cbegin(), planIdcs.cend(),
				   aggregateExchangeMode, leafLevel);
	  }

      }

      friend std::ostream& operator<< (std::ostream &out, const Plan &e)
      {
	out << "(level: " << e.level
	    << ", nNodes " << e.planIdcs.size()
	    << ", aggregateExchangeMode " << e.aggregateExchangeMode << ")";
	return out;
      }

      auto getLevel() const
      {
	return level;
      }

      auto getAggregateExchangeMode() const
      {
	return aggregateExchangeMode;
      }
      
    private:
      const uint32_t level;
      const bool aggregateExchangeMode;
      std::vector<uint32_t> planIdcs;
      // required for specific mode
      std::vector<std::tuple<std::vector<IDX>, std::vector<IDX>, size_t>>
      inBoxes_planProbList;
    };
    
    class Plans
    {
      static const bool cleanUpAfterHighLevelExchange=false;
      using IDX=uint32_t;
      using DIM=V2<IDX>;
    public:      
      
      template<typename QT>
      static auto regularTileIndices(const size_t nTilesFull)
      {
	// using IDX=uint32_t;
	// using DIM=V2<IDX>;		

	const auto exp4Upper = helper::ilog4_ceil(nTilesFull);
	const auto exp4Lower = helper::ilog4_floor(nTilesFull);
	
	const auto nLeavesUpper=
	  helper::ipow4(exp4Upper);

	const auto nLeavesLower=
	  helper::ipow4(exp4Lower);

	const DIM gridDimLower(helper::ipow2(exp4Lower),
			       helper::ipow2(exp4Lower));

	const DIM gridDimUpper(helper::ipow2(exp4Upper),
			       helper::ipow2(exp4Upper));

	assert(helper::ii2n(gridDimLower)==nLeavesLower);
	assert(helper::ii2n(gridDimUpper)==nLeavesUpper);

	hassertm2(nLeavesLower <= nTilesFull, nLeavesLower, nTilesFull);
	hassertm2(nLeavesUpper >= nTilesFull, nLeavesUpper, nTilesFull);

	

	const auto grid2qtUpper = genMap_grid2qt(gridDimUpper);
	
	if(nLeavesUpper==nLeavesLower)
	  return std::make_tuple(helper::range_n<std::vector<IDX>>(nLeavesUpper),grid2qtUpper);

	QT qt(nLeavesUpper);
	
	// initialize with upper left quadrant
	auto idcs=helper::range_n<std::vector<IDX>>(nLeavesLower);
	idcs.reserve(nTilesFull);

	auto addIdx = [&idcs](auto idx)
	{
	  hassertm(std::find(idcs.begin(), idcs.end(), idx)==idcs.end(),
		   idx);
	  idcs.push_back(idx);
	};
	
	
	//
	// extend to top-right quadrant
	//
	{
	  const auto offIdcs = helper::range2_n(gridDimLower);
	  for(size_t i=0;
	      i<std::min(offIdcs.size(), nTilesFull-nLeavesLower);
	      i++)
	    {
	      // swap x and y here!!
	      const auto off=DIM(offIdcs[i].y, offIdcs[i].x);
	      const auto idx=grid2qtUpper
		[helper::ii2i(DIM(gridDimLower.x, 0)+off, gridDimUpper)];
	      
	      addIdx(idx);
	    }
	}

	//
	// extend to bottom half
	//
	if(idcs.size()<nTilesFull)
	  {
	    auto offIdcs = helper::range2_n
	      (DIM(gridDimUpper.x, gridDimLower.y));	  

	    assert(idcs.size()==2*nLeavesLower);
	    const auto nIdcsRemaining = nTilesFull-2*nLeavesLower;

	    assert(nIdcsRemaining < offIdcs.size());

	    offIdcs.resize(nIdcsRemaining);
	  
	    for(const auto & off : offIdcs)
	      {	      
		const auto idx=grid2qtUpper
		  [helper::ii2i(DIM(0, gridDimLower.y)+off, gridDimUpper)];
		addIdx(idx);
	      }
	  }	
	
	assert(idcs.size()==nTilesFull);
	return std::make_tuple(idcs, grid2qtUpper);
      }

      template<typename QT>
      static auto isRegularMap(const size_t nTilesFull, const size_t nTiles)
      {
	//QT qt(nTilesFull);
	std::vector<bool> isRegular(nTiles, false);
	{
	  auto regularLeafIdcs =
	    std::get<0>(regularTileIndices<QT>(nTilesFull));

	  assert(regularLeafIdcs.size()==nTilesFull);
	      
	      
	  for(const auto & idx : regularLeafIdcs)
	    {
	      assert(idx<isRegular.size());
	      assert(!isRegular[idx]);
	      isRegular[idx]=true;
	    }
	}
	return isRegular;
      }
      
      template<typename QT>
      Plans(const QT& qt, size_t nNodesLevel, size_t chunkSizeMax_in, size_t maxNodeExchangeLevel) :
	chunkSizeMax(chunkSizeMax_in)
      {
	using IDX=uint32_t;
	auto regularNodeIndices = [&]
	  (const std::vector<IDX>& regularLeafIdcs,
	   const uint32_t level)
	{
	  const auto exp4Upper
	    = helper::ilog4_ceil(regularLeafIdcs.size());	

	  QT qt(helper::ipow4(exp4Upper));

	  std::vector<bool> isRegular(qt.nLeaves(), false);
	  for(const auto & idx : regularLeafIdcs)
	    {
	      assert(idx<isRegular.size());
	      isRegular[idx]=true;
	    }
	
	
	  const auto qt2leavesUpper = genMap_qt2leaves(qt);

	  // check if all leaves are regular
	  const auto levelOffset = qt.getLevelOffset(level);
	  const auto nNodesLevel = qt.nNodesLevel(level);

	  std::vector<IDX> idcs;
	
	  for(const auto & levelIdxNode : helper::range_n(nNodesLevel))
	    {
	      const auto & leaves = qt2leavesUpper[levelOffset+levelIdxNode];

	      bool allRegular=true;
	      for(const auto & leaf : leaves)
		{
		  assert(leaf < isRegular.size());
		  allRegular = allRegular && isRegular[leaf];
		}
	      if(allRegular)
		idcs.push_back(levelIdxNode);	    
	    }
	  return idcs;
	};

	// std::vector<IDX> regularLeafIdcs;
	// std::vector<DIM> grid2qt;
	// QT qt;
	const auto [regularLeafIdcs, grid2qt]
	  = regularTileIndices<QT>(nNodesLevel);
	
	
	size_t level = 0;
		
	//size_t nLeavesLevel=1;

	//const auto nLevels = qt.nLevels(nNodesLevel/N);
	const auto nLevels = qt.nLevels_ceil(nNodesLevel);
	
	//while(nNodesLevel>=N)
	//while(nNodesLevel>chunkSizeMax*32)
	//while(nNodesLevel>chunkSizeMax/**4*/)
	while(nNodesLevel>chunkSizeMax && maxNodeExchangeLevel>=level)
	  {
	    auto addPlan = [&](const auto& plan)
	    {	      
	      assert(level<nLevels);

	      const size_t from_n = plansProbs.empty() ? 0 : plansProbs.back().second;

	      const size_t shift = /*2**/(nLevels-1-level);
	      assert(shift < 64);
#if 0
	      const size_t n = shift;
#else
	      const size_t n = 1<<shift;
#endif
	      std::cout << "level "<< level << " vs nLevels "  << nLevels << " | n: " << n << " " << from_n+n << std::endl;

	      plansProbs.emplace_back(plan, from_n+n);
	    };
	      
	    if(level==0)
	      {		
		Plan plan(level, false, regularLeafIdcs, grid2qt);
		addPlan(plan);
		plansReorderChildren.push_back(regularLeafIdcs);
	      }
	    if(level>0)
	      {
		const auto rni = regularNodeIndices(regularLeafIdcs, level);
		Plan plan(level, true, rni, grid2qt);
		addPlan(plan);
		plansReorderChildren.push_back(rni);
	      }

	    // std::cout << "plansReorderChildren " << level << ": \n";
	    // for(const auto e : plansReorderChildren.back())
	    //   std::cout << e << " ";
	    // std::cout << std::endl;
	      
	    level++;
	    nNodesLevel/=4;
	  }
	std::cout << "there are " << plansProbs.size() << " plans for different levels (maxNodeExchangeLevel: "<< maxNodeExchangeLevel<< ")\n";
	assert(!plansProbs.empty());
      }
      
      template<typename RNG>
      auto operator()(RNG& rng)
      {
	std::vector<uint32_t>::const_iterator plan_begin, plan_end;
	bool aggregateExchangeMode=false;
	size_t leafLevel=-1;
	size_t chunkSize=-1;

	
	if(prevLevel==0 || !cleanUpAfterHighLevelExchange)
	  {
	    
	    std::uniform_int_distribution<size_t> dis(static_cast<size_t>(0), plansProbs.back().second);

	    const size_t n=dis(rng);

	    for(auto & p : plansProbs)
	      //#warning "MODIFICATION FOR TESTING EXCHANGES WITH LEVEL=5"
	      if(n<=p.second/*p.first.getLevel()==5*/)
		{		  
		  std::tie(plan_begin, plan_end, aggregateExchangeMode, leafLevel)
		    = p.first(rng);
		  break;
		}
	    
	    
	    
	    prevLevel=leafLevel;
	    chunkSize=chunkSizeMax;
	  }
	else
	  {
	    assert(plansProbs.front().first.getLevel()==0);
	    prevLevel--;

	    leafLevel=prevLevel;
	    aggregateExchangeMode=(leafLevel>0);
	    

	    // std::tie(plan, aggregateExchangeMode, leafLevel)	 
	    //   = plansProbs[prevLevel].first(rng);

	    // assert(prevLevel==leafLevel);

	    

	    chunkSize=4;
	    //plan=helper::range_n<std::vector<uint32_t>>(plan.size());
	    plan_begin=plansReorderChildren[prevLevel].cbegin();
	    plan_end=plansReorderChildren[prevLevel].cend();
	
	  }
	
	return std::make_tuple(plan_begin, plan_end, aggregateExchangeMode, leafLevel, chunkSize);
      }
      
    private:

      std::vector<std::pair<Plan, size_t>> plansProbs;
      size_t chunkSizeMax;
      size_t prevLevel=0;

      std::vector<std::vector<uint32_t>> plansReorderChildren;
    };
                
  }
}

#endif //__SUPERTILES_PLACE_PLAN__
