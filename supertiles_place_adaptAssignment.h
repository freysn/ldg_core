#ifndef __SUPERTILES_PLACE_ADAPT_ASSIGNMENT__
#define __SUPERTILES_PLACE_ADAPT_ASSIGNMENT__

namespace supertiles
{
  namespace place
  {

    template<typename QT_IT, typename NodeLeafCounts, typename D, typename D_TILE>
    void addRm(QT_IT qt_it_add, QT_IT qt_it_rm, D* supertiles, NodeLeafCounts& nodeLeafCounts,
	       const D_TILE* tile, size_t nElemsTile)
    {
      while(qt_it_add.idx != qt_it_rm.idx)
	{
	  for(const auto & i : helper::range_n(nElemsTile))
	    {
	      
	      const auto & v = tile[i];
	      supertiles[qt_it_add.idx*nElemsTile+i]+= v;

	      
	      
	      auto & e_rm = supertiles[qt_it_rm.idx*nElemsTile+i];
	      e_rm-= v;	      

	      //hassertm2(e_rm >= -m_eps, e_rm, qt_it_add.nElemsLevel);
	    }

	  nodeLeafCounts.addNode(qt_it_add.idx);
	  nodeLeafCounts.rmNode(qt_it_rm.idx);
	  
	  {
	    const bool success0 = qt_it_add();
	    const bool success1 = qt_it_rm();

	    hassert(success0==success1);
	  }
	}
    }

    template<typename IT, typename N, typename QT, typename A, typename M>
    auto adaptAggregate2Leaf(IT idcs_begin, IT idcs_end,  const N& nodes2leaves, const QT& qt, A& qtLeafAssignment_out, /*const A& qtNodeAssignment_vanilla,*/ const M& movedFromQTPos, size_t supertileLevelOffset)
    {

	//
	// convert higher level node assigment to leaf assignments
	//
      //auto qtLeafAssignment_out = qtLeafAssignment;
      const auto qtLeafAssignment = qtLeafAssignment_out;
      
      auto movedFromQTPos_out = helper::range_n<std::vector<uint32_t>>(qt.nLeaves());

	std::vector<uint32_t> leafIdcs;
	if(idcs_begin != idcs_end)
	  {
	    assert(nodes2leaves[supertileLevelOffset+(*idcs_begin)].size()
		   ==
		   nodes2leaves[supertileLevelOffset].size());
	    
	    const auto n = nodes2leaves[supertileLevelOffset+(*idcs_begin)].size()*(idcs_end-idcs_begin);
	    leafIdcs.reserve(n);
	  }
		      
	//for(size_t qtPosIdx=0; qtPosIdx<qtNodeAssignment_vanilla.size(); qtPosIdx++)
	for(auto it=idcs_begin; it != idcs_end; it++)
	  {	    
	    const auto qtPosIdx=(*it);
	    const auto nodeId =qtPosIdx+supertileLevelOffset;

	    
	    const auto oldQTPosIdx = movedFromQTPos[qtPosIdx];
	    
	    const auto oldNodeId = oldQTPosIdx+supertileLevelOffset;
			  
	    //std::cout << "[" << qtPosIdx << "|" << nodeId << "|" << oldQTPosIdx<< "]\n";

	    //hassertm4(qtNodeAssignment_vanilla[qtPosIdx]==movedFromQTPos[qtPosIdx], qtPosIdx, qt.nNodesLevel(5), qtNodeAssignment_vanilla[qtPosIdx],movedFromQTPos[qtPosIdx]);

	    if(qtPosIdx != oldQTPosIdx/* || true*/)
	      {
		const auto & leavesTo = nodes2leaves[nodeId];
		const auto & leavesFrom = nodes2leaves[oldNodeId];

		assert(leavesTo.size()==leavesFrom.size());
			      
		for(const auto & leafId : helper::range_n(leavesTo.size()))
		  {
		    const auto leafTo = leavesTo[leafId];
		    const auto leafFrom = leavesFrom[leafId];

		    assert(leafFrom != leafTo);
		    const auto tileId = qtLeafAssignment[leafFrom];
				  
		    qtLeafAssignment_out[leafTo] = tileId;
		    movedFromQTPos_out[leafTo] = leafFrom;

		    leafIdcs.push_back(leafTo);
		  }
	      }
	  }
	//std::cout << std::endl;

	//qtLeafAssignment = qtLeafAssignment_out;

// #ifndef NDEBUG
	// THIS ONLY HOLDS TRUE IF NOT FILTERING qtPosIdx != oldQTPosIdx
// 	RangeLeaves<IT> rangeLeaves(idcs_begin, idcs_end, supertileLevelOffset, nodes2leaves);
// 	hassertm2(leafIdcs.size()==rangeLeaves.size(), leafIdcs.size(), rangeLeaves.size());
// 	for(size_t i=0;i<leafIdcs.size();i++)
// 	  hassertm4(leafIdcs[i]==rangeLeaves.size(), leafIdcs[i], rangeLeaves[i], i, leafIdcs.size());
// #endif

	return std::make_tuple(leafIdcs, movedFromQTPos_out);
    }


    template<typename IT, typename S, typename T, typename NodeLeafCounts, typename QT, typename A, typename M>
    void adaptAssignmentLeaf(IT idcs_begin, IT idcs_end, S& supertiles, const T& tiles,
			     NodeLeafCounts& nodeLeafCounts,
			     const QT& qt, const size_t nElemsTile,  A& qtLeafAssignment, M& movedFromQTPos)
    {
      //for(size_t qtPosIdx=0; qtPosIdx<qtLeafAssignment.size(); qtPosIdx++)
      //std::cout << "---------\n";
      for(auto it=idcs_begin; it != idcs_end; it++)
	{
	  const auto qtPosIdx=*it;
	  const auto tileId = qtLeafAssignment[qtPosIdx];
	  const auto oldQTPosIdx = movedFromQTPos[qtPosIdx];
	  const bool isVoid = (tileId==voidTileIdx);
	  
	  //std::cout << "qtPosIdx " << qtPosIdx << " tileId " << tileId << " void: " << isVoid << " nodeLeafCounts " << nodeLeafCounts(qtPosIdx) << std::endl;

	  if(qtPosIdx != oldQTPosIdx && !isVoid)
	    {		  
	      addRm(qt(qtPosIdx), qt(oldQTPosIdx),
		    &supertiles[0], nodeLeafCounts,
		    /*isVoid ? static_cast<D*>(0) :*/ &tiles[nElemsTile*tileId],
		    nElemsTile);
	    }
	}
    }

    template<typename IT, typename S, typename T, typename N, typename NodeLeafCounts, typename QT, typename A, typename M>
    void adaptAssignment(const IT idcs_begin, const IT idcs_end, S& supertiles, const T& tiles, const N& nodes2leaves,
			 NodeLeafCounts& nodeLeafCounts, const QT& qt, const size_t nElemsTile,  A& qtLeafAssignment,  const M& movedFromQTPos, bool aggregateExchangeMode, size_t supertileLevelOffset)
    {
#ifndef NDEBUG
      nodeLeafCounts.dbg_sanityCheck();
#endif

      auto _adaptAssignmentLeaf = [&supertiles, &tiles, &nodeLeafCounts, &qt, &nElemsTile, &qtLeafAssignment]
	(const auto idcs_begin, const auto idcs_end, const auto & movedFromQTPos)
      {
	adaptAssignmentLeaf(idcs_begin, idcs_end, supertiles, tiles, nodeLeafCounts, qt, nElemsTile,  qtLeafAssignment,
			    movedFromQTPos);
      };
      
      if(aggregateExchangeMode)
	{
	  const auto [_idcs, _movedFromQTPos]
	    = adaptAggregate2Leaf(idcs_begin, idcs_end, nodes2leaves, qt, qtLeafAssignment, movedFromQTPos, supertileLevelOffset);
	  //adaptAssignmentLeaf(_idcs.begin(), _idcs.end(), supertiles, tiles, nodeLeafCounts, qt, nElemsTile,  qtLeafAssignment, _movedFromQTPos);
	  _adaptAssignmentLeaf(_idcs.begin(), _idcs.end(), _movedFromQTPos);
	  
	}
      else
	//adaptAssignmentLeaf(idcs_begin, idcs_end, supertiles, tiles, nodeLeafCounts, qt, nElemsTile,  qtLeafAssignment, movedFromQTPos);
	_adaptAssignmentLeaf(idcs_begin, idcs_end, movedFromQTPos);

#ifndef NDEBUG
      nodeLeafCounts.dbg_sanityCheck();
#endif
    }

#if 0
    template<typename IT, typename S, typename T, typename N, typename QT, typename A, typename M>
    void DEPRECATED_adaptAssignment(IT idcs_begin, IT idcs_end, S& supertiles, const T& tiles, const N& nodes2leaves, const QT& qt, const size_t nElemsTile,  A& qtLeafAssignment, A& qtNodeAssignment_vanilla, M& movedFromQTPos, bool aggregateExchangeMode)
    {
      if(aggregateExchangeMode)
	{
	  //
	  // convert higher level node assigment to leaf assignments
	  //
	  auto qtLeafAssignment_out = qtLeafAssignment;
	  auto movedFromQTPos_out = helper::range_n(qt.nLeaves());

		      
		      
	  //for(size_t qtPosIdx=0; qtPosIdx<qtNodeAssignment_vanilla.size(); qtPosIdx++)
	  for(auto it=idcs_begin; it != idcs_end; it++)
	    {
	      const auto qtPosIdx=*it;
	      const auto nodeId =qtPosIdx;
			  
	      const auto oldQTPosIdx = movedFromQTPos[qtPosIdx];
			  
	      //std::cout << "[" << qtPosIdx << "|" << nodeId << "|" << oldQTPosIdx<< "]\n";

	      assert(qtNodeAssignment_vanilla[qtPosIdx]==movedFromQTPos[qtPosIdx]);
			  
	      if(qtPosIdx != oldQTPosIdx)
		{
		  const auto & leavesTo = nodes2leaves[nodeId];
		  const auto & leavesFrom = nodes2leaves[oldQTPosIdx];

		  assert(leavesTo.size()==leavesFrom.size());
			      
		  for(const auto & leafId : helper::range_n(leavesTo.size()))
		    {
		      const auto leafTo = leavesTo[leafId];
		      const auto leafFrom = leavesFrom[leafId];

		      assert(leafFrom != leafTo);
		      const auto tileId = qtLeafAssignment[leafFrom];
				  
		      qtLeafAssignment_out[leafTo] = tileId;
		      movedFromQTPos_out[leafTo] = leafFrom;
		    }
		}
	    }
	  //std::cout << std::endl;

	  qtLeafAssignment = qtLeafAssignment_out;
	  movedFromQTPos = movedFromQTPos_out;
	}

      for(size_t qtPosIdx=0; qtPosIdx<qtLeafAssignment.size(); qtPosIdx++)
	{
	  const auto tileId = qtLeafAssignment[qtPosIdx];
	  const auto oldQTPosIdx = movedFromQTPos[qtPosIdx];
	  if(qtPosIdx != oldQTPosIdx)
	    {		  
	      addRm(qt(qtPosIdx), qt(oldQTPosIdx),
		    &supertiles[0], &tiles[nElemsTile*tileId], nElemsTile);
	    }
	}	      		    
    }
#endif
  }
}

#endif //__SUPERTILES_PLACE_ADAPT_ASSIGNMENT__
