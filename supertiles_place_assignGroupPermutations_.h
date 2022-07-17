
template<distFuncType_t DIST_FUNC,  typename ASS, typename D_IT, typename IDX, typename NodeLeafCounts, typename NodeLeafCounts0, typename N2L, typename QT, typename P>
auto assignGroupPermutations
(ASS& qtLeafAssignment,
 std::vector<size_t>& movedFromQTPos,
 const D_IT supertiles,
 const D_IT supertiles_0,
 //const std::vector<D>& tiles,
 const size_t nElemsTile,
 const NodeLeafCounts& nodeLeafCounts,
 const NodeLeafCounts0& nodeLeafCounts_0,
 const N2L& nodes2leaves_0,
 const QT& qt,
 const std::vector<IDX>& qtIdxv, const std::vector<uint32_t>& qtNeighbors,
 const size_t nNeighbors, const double neighborFac,
 const P& assignmentCandidatesv)
{

#ifndef __NO_OMP
  // currnt cost cachign only works when no openmp is enabled
  assert(false); 
#endif

  using node_t=uint32_t;
  
  using D=double;

  
  //const size_t maxNLevelsUp=-1;
  const size_t maxNLevelsUp=1;
  const bool beginOneLevelUp=true;
  
  const auto N = qtIdxv.size();

  assert(N>0);
  if(N==1)
    return false;

  hassertm2(N < assignmentCandidatesv.size(), N, assignmentCandidatesv.size());
  const auto & assignmentCandidates=assignmentCandidatesv[N];
  assert(assignmentCandidates.size()==factorial(N));   
   
   
  using QT_IT = supertiles_QuadTree_iterator_up<IDX>;
   
  std::vector<IDX> tileIdxv(N);
  std::vector<QT_IT> qt_itv(N);

   
   
   
  for(uint32_t i=0; i<N; i++)
    {       
      const auto & p = qtIdxv[i];
      hassertm2(p<qtLeafAssignment.size(), p, qtLeafAssignment.size());
      tileIdxv[i] = qtLeafAssignment[p];
      qt_itv[i] = qt(p);
      assert(qtLeafAssignment[qtIdxv[i]] == tileIdxv[i]);
      assert(movedFromQTPos[p]==p);	 
    }
  //std::cout << std::endl;

  const size_t nThreads =
#ifndef __NO_OMP
    omp_get_max_threads();
#else
  1;
#endif

  using idcs_t = typename
    std::remove_const<typename std::remove_reference<decltype(assignmentCandidates[0])>::type>::type;
   
  std::vector<std::pair<D, idcs_t>>
    bestv(nThreads, std::make_pair(std::numeric_limits<D>::max(), idcs_t()));

  D minCost = std::numeric_limits<D>::max();
  
#ifndef __NO_OMP
#pragma omp parallel for
#endif
  for(size_t assignmentIdx=0; assignmentIdx<assignmentCandidates.size();
      assignmentIdx++)
    {
      const auto threadIdx =
#ifndef __NO_OMP
	omp_get_thread_num();
#else
      0;
#endif

      D minCost_;
#ifndef __NO_OMP
      #pragma omp atomic read
#endif
      minCost_ = minCost;

      auto & best = bestv[threadIdx];
      
#ifndef __NO_OMP
      minCost_ = std::min(minCost_, best.first);
#else
      hassertm2(minCost_==best.first, minCost_, best.first);
#endif
      
      
       
      const auto & assignment = assignmentCandidates[assignmentIdx];
      std::vector<size_t> replaceLeaves(N);

      for(const auto i : helper::range_n(N))
	replaceLeaves[assignment[i]] = qtIdxv[i];
     
      D new_cost=0.;

      auto stillBetter = [&new_cost, &minCost_, &best]()
      {return new_cost < std::min(minCost_, best.first);};

      
      for(uint32_t qtPosListIdx=0; qtPosListIdx<N && stillBetter();
	  qtPosListIdx++)
	{	      
	  auto qt_it_ref = qt_itv[qtPosListIdx];

	  CostVecCache costCache;
	  new_cost
	    += costVec<costVecMode_replace, DIST_FUNC>
	    (costCache,
	     qt_itv, qt_it_ref, nodes2leaves_0[replaceLeaves[qtPosListIdx]],
	     &supertiles[0],
	     &supertiles_0[0],
	     nElemsTile, nodeLeafCounts, nodeLeafCounts_0, maxNLevelsUp, beginOneLevelUp, replaceLeaves);

	  for(uint32_t i=0; i<nNeighbors
		&& (new_cost < best.first) && stillBetter();
	      i++)
	    {
	      const auto neighbor = qtNeighbors[nNeighbors*qt_it_ref.idx+i];
	      
	      if(neighbor != -1)
		{
		  new_cost +=
		    neighborFac*
		    costVec<costVecMode_replace, DIST_FUNC>
		    (costCache,
		     qt_itv, qt(neighbor),
		     nodes2leaves_0[replaceLeaves[qtPosListIdx]],
		     &supertiles[0],
		     &supertiles_0[0],
		     nElemsTile, nodeLeafCounts, nodeLeafCounts_0, maxNLevelsUp,
		     beginOneLevelUp, replaceLeaves);
		}
	    }	  
	}
       
      if(stillBetter())
	{
	  best=std::make_pair(new_cost, assignment);
	

      // it may happen in parallel execution
      // that minCost_ > best.first

      // NOTE THAT THIS IS ONLY FOR EARLY TERMINATION
      // THIS DOES NOT AFFECT THE OUTCOME (THIS IS COLLECTED BELOW)
	 
	  D minCost__;
	  
#ifndef __NO_OMP
#pragma omp atomic read
#endif
	  minCost__ = minCost;

	  if(new_cost < minCost__)
	    {
#ifndef __NO_OMP
#pragma omp atomic write
#endif
	      minCost = new_cost;

	    }	
	}
    }

  auto & best = bestv[0];
   
  for(size_t threadIdx=1; threadIdx<nThreads; threadIdx++)
    if(bestv[threadIdx].first < best.first)
      best = bestv[threadIdx];
   
  const bool adapt = !std::is_sorted(std::get<1>(best).begin(),
				     std::get<1>(best).end());
   
  if(adapt)
    {
      auto assignment = std::get<1>(best);
      // auto updateAssignment = [&assignment, &qtIdxv, &tileIdxv, &N]
      // 	 (auto& qtLeafAssignment, auto& movedFromQTPos)
      {
	for(const auto i : helper::range_n(N))
	  {
	    const auto newQTPos = qtIdxv[assignment[i]];
	    const auto oldQTPos = qtIdxv[i];

	    assert(newQTPos < qtLeafAssignment.size());
	    qtLeafAssignment[newQTPos] = tileIdxv[i];
	    if(!movedFromQTPos.empty())
	      {
		assert(newQTPos < movedFromQTPos.size());
		movedFromQTPos[newQTPos] = oldQTPos;
	      }
	  }
      };
       
      //updateAssignment(qtLeafAssignment, movedFromQTPos);
    }
  return adapt;

}
