
template<distFuncType_t DIST_FUNC, typename MIDX, typename ASS, typename D_IT, typename IDX, typename NodeLeafCounts, typename NodeLeafCounts0, typename N2L, typename QT/*, typename P*/>
auto assignGroupPermutations
(ASS& qtLeafAssignment,
 std::vector<MIDX>& movedFromQTPos,
 const D_IT supertiles,
 const D_IT supertiles_0,
 //const std::vector<D>& tiles,
 const size_t nElemsTile,
 const NodeLeafCounts& nodeLeafCounts,
 const NodeLeafCounts0& nodeLeafCounts_0,
 const N2L& nodes2leaves_0,
 const QT& qt,
 const std::vector<IDX>& qtIdxv, const std::vector<uint32_t>& qtNeighbors,
 const size_t nNeighbors, const double neighborFac/*,
						    const P& assignmentCandidatesv*/, const size_t maxNLevelsUp)
{  
  
  using D=double;

  
  //const size_t maxNLevelsUp=-1;
  //const size_t maxNLevelsUp=1;
  //const bool beginOneLevelUp=true;
  
  const auto N = qtIdxv.size();

  assert(N>0);
  if(N==1)
    return false;
   
   
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
  // using idcs_t = typename
  //   std::remove_const<typename std::remove_reference<decltype(assignmentCandidates[0])>::type>::type;
   
  // std::vector<std::pair<D, idcs_t>>
  //   bestv(nThreads, std::make_pair(std::numeric_limits<D>::max(), idcs_t()));


  std::vector<D> costMatrix(N*N, 0.);

  std::vector<CostVecCache<D>> costCachev(N);
  
  for(size_t x=0; x<costMatrix.size(); x++)
    {      
      const auto qtPosListIdx = x/N;
      const auto refNodeIdx = x-qtPosListIdx*N;

      auto & cost = costMatrix[helper::hungarian_costMatrixIdx(refNodeIdx,qtPosListIdx,N)];
      
      const auto & refLeafNodeIdcs=nodes2leaves_0[qtIdxv[refNodeIdx]];
      
      auto qt_it_ref = qt_itv[qtPosListIdx];

      auto & costCache=costCachev[refNodeIdx];

      // can begin one level up as base level is element itself by definition

      // std::cout << "test: " <<
      // 	costVec<costVecMode_replace, DIST_FUNC>
      // 	(costCache,
      // 	 qt_itv, qt_it_ref,
      // 	 refLeafNodeIdcs,
      // 	 &supertiles[0],
      // 	 &supertiles_0[0],
      // 	 nElemsTile, nodeLeafCounts, nodeLeafCounts_0, 0, 0)
      // 		<< std::endl;
      
	cost
	+= costVec<costVecMode_replace, DIST_FUNC>
	(costCache,
	 qt_itv, qt_it_ref,
	 refLeafNodeIdcs,
	 &supertiles[0],
	 &supertiles_0[0],
	 nElemsTile, nodeLeafCounts, nodeLeafCounts_0, maxNLevelsUp, 1);

      for(uint32_t i=0; i<nNeighbors;i++)
	{
	  const auto neighbor = qtNeighbors[nNeighbors*qt_it_ref.idx+i];
	      
	  if(neighbor != -1)
	    {
	      cost
		+=
		neighborFac*
		costVec<costVecMode_replace, DIST_FUNC>
		(costCache,
		 qt_itv, qt(neighbor),
		 refLeafNodeIdcs,
		 &supertiles[0],
		 &supertiles_0[0],
		 nElemsTile, nodeLeafCounts, nodeLeafCounts_0, maxNLevelsUp,
		 0);
	    }
	}
    }

  std::vector<IDX> assignment(N, -1);
  helper::hungarian_assignmentoptimal
    (&assignment[0], &costMatrix[0], static_cast<IDX>(N), static_cast<IDX>(N));
   
  const bool adapt = !std::is_sorted(assignment.begin(), assignment.end());
   
  if(adapt)
    {
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
