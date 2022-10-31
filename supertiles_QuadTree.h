#ifndef __SUPERTILES_QUAD_TREE__
#define __SUPERTILES_QUAD_TREE__

#include "helper_util.h"
#include "helper_assert.h"
#include "helper_imath.h"
#include "supertiles_place_NodeLeafCounts.h"

namespace supertiles
{
namespace place
{
  constexpr uint32_t voidTileIdx=-1;
}
};

template<typename IDX>
class supertiles_QuadTree_iterator_up
{
  static const IDX idx_void = supertiles::place::voidTileIdx;;
public:
  supertiles_QuadTree_iterator_up()
  {}
    
  supertiles_QuadTree_iterator_up(IDX idx_in,
				  IDX nElemsLevel_in)
  {
    hassertm3(idx_in==uint32_t(-1) || idx_in < nElemsLevel_in, idx_in, nElemsLevel_in, uint32_t(-1));

    if(idx_in==uint32_t(-1))
      idx = idx_void;
    else
      idx = idx_in;
    nElemsLevel=nElemsLevel_in;
    levelOffset=0;
  }
    
  bool operator()()
  {
    assert(nElemsLevel > 0);
    auto levelOffsetNext = levelOffset+nElemsLevel;
    nElemsLevel/=4;
    if(nElemsLevel==0)
      return false;

    if(idx != idx_void)
      idx = (idx-levelOffset)/4 + levelOffsetNext;
    levelOffset = levelOffsetNext;
    return true;
  }

  /*
  IDX getChild0Idx() const
  {
    assert(idx >= levelOffset);
    auto levelIdx = idx-levelOffset;

    levelIdx *= 4;

    assert(levelOffset>=4*nElemsLevel);
    const auto levelOffsetPrev = levelOffset-4*nElemsLevel;

    return levelOffsetPrev+levelIdx;    
  }
  */
  
  bool isVoid() const
  {
    return idx==idx_void;
  }
  
  IDX idx=-1;
  IDX nElemsLevel=-1;
  IDX levelOffset=-1;
};

template<typename T>
inline std::ostream& operator<< (std::ostream &out, const supertiles_QuadTree_iterator_up<T> &e)
{
  out << "(idx: " << e.idx << ", nElemsLevel " << e.nElemsLevel << ", levelOffset " << e.levelOffset << ")";
  return out;
}


template<typename IDX>
class supertiles_QuadTree
{
public:  

  supertiles_QuadTree()
  {
    nElemsLevel_begin = -1;
  }

  static auto nLeavesPerNode(size_t level)
  {
    return 1<<(2*level);
  }
  
  auto nNodesLevel(size_t level) const
  {
    return nElemsLevel_begin/nLeavesPerNode(level);
  }

  supertiles_QuadTree<IDX> higherLevelCopy(size_t level) const
  {
    hassertm2(level < nLevels(), level, nLevels());
    return supertiles_QuadTree(nNodesLevel(level));
  }

  
  auto getLevelOffset(size_t level) const
  {
    auto it = (*this)(0);
    for(size_t l=0; l<level; l++)
      {
	const bool success = it();
	hassert(success);
      }
    return it.levelOffset;
  }

  static IDX nLeavesMax(IDX nElems)
  {
    IDX level = std::ceil((log(nElems) / log(4)));
    return std::pow(4, level);
  }
  
  supertiles_QuadTree(IDX nElems)
  {
    assert(nElems>0);
    nElemsLevel_begin = nLeavesMax(nElems);
    //std::cout << "generate quadtree with " << nElemsLevel_begin << " leaves for " << nElems << " tiles\n";
  }

  template<typename I>
  supertiles_QuadTree(const V2<I> gridDim)
  {    
    nElemsLevel_begin = helper::ii2n(gridDim);
  }

  IDX firstChild(IDX idx) const
  {
    if(idx < nElemsLevel_begin)
      return -1;

    const auto level=getLevel(idx);

    const auto o=(idx-getLevelOffset(level))*4+getLevelOffset(level-1);
    assert(parent(o)==idx);
    return o;
  }
  
  supertiles_QuadTree_iterator_up<IDX> operator()(IDX idx) const
  {
    //assert(idx<nElemsLevel_begin);
    

    if(idx<nElemsLevel_begin)
      return supertiles_QuadTree_iterator_up<IDX>(idx, nElemsLevel_begin);

       supertiles_QuadTree_iterator_up<IDX> qt_it(0, nElemsLevel_begin);

    const auto level=getLevel(idx);
    
    for(IDX l=0; l<level; l++)
      qt_it();

    assert(idx>=qt_it.idx);
    qt_it.idx=idx;
    
    return qt_it;
  }

  IDX parent(IDX idx) const
  {
    auto qt_it=(*this)(idx);
    if(!qt_it())
      return -1;
    return qt_it.idx;
  }

  IDX nElems() const
  {
    //std::cout << "---------------------\n";
    auto it = supertiles_QuadTree_iterator_up<IDX>(0, nElemsLevel_begin);
    //std::cout << it << std::endl;
    while(it())
      {
	//std::cout << it << std::endl;
      }
    //std::cout << "###############\n";
    return it.levelOffset+1;
  }

  IDX nElemsLevel(uint32_t level) const
  {
    //std::cout << "---------------------\n";
    auto it = supertiles_QuadTree_iterator_up<IDX>(0, nElemsLevel_begin);
    //std::cout << it << std::endl;
    
    for(IDX l=0; l<level; l++)
      {
	if(!it())
	  return 0;
	//std::cout << it << std::endl;
      }
    //std::cout << "###############\n";
    return it.nElemsLevel;
  }

  IDX nNodes() const
  {
    return nElems();
  }

  IDX nLeaves() const
  {
    return nElemsLevel_begin;
  }

  IDX nLevels() const
  {
    return nLevels(nLeaves());
    // const auto v = log(nLeaves()) / log(4);    
    // return v + 1. + .5;
  }

  static IDX nLevels(size_t nLeaves)
  {
    assert(helper::isPow4(nLeaves));

    const IDX rslt = nLeaves==1 ? 1 : helper::ilog4_ceil(nLeaves)+1;

#ifndef NDEBUG
    const auto v = log(nLeaves) / log(4); 
    const IDX rslt_test = v + 1. + .5;    
    hassertm3(rslt==rslt_test, rslt, rslt_test, nLeaves);
#endif
    
    return rslt;
  }

  static IDX nLevels_floor(size_t nLeaves)
  {
    return helper::ilog4_floor(nLeaves)+1;
  }

  static IDX nLevels_ceil(size_t nLeaves)
  {
    return helper::ilog4_ceil(nLeaves)+1;
  }


  IDX getLevel(IDX node) const
  {
    assert(node!=static_cast<IDX>(-1));
    hassertm2(node < nElems(), node, nElems());

    IDX nElemsLevel = nLeaves();
    IDX level=0;
    while(node>=nElemsLevel)
      {
	assert(nElemsLevel > 0);
	node -=nElemsLevel;
	nElemsLevel /= 4;
	level++;
      }

    assert(level < nLevels());
    return level;
  }

  static IDX end_idx()
  {
    return static_cast<IDX>(-1);
  }
  
  IDX nElemsLevel_begin;  
};


template<typename DIM>
auto genMap_supergrid2qt(DIM inDim)
{
  using pxR_t = std::pair<DIM, DIM>;

  using idxR_t = std::pair<uint32_t, uint32_t>;

  std::vector<std::pair<pxR_t, idxR_t>>
    v(1, std::make_pair(pxR_t(DIM(0,0), inDim), idxR_t(0, helper::ii2n(inDim))));

  auto v_levels = v;
    
  while(v.front().second.second>1)
    {
      std::vector<std::pair<pxR_t, idxR_t>> vo;
      std::swap(v, vo);
      for(const auto & e : vo)
	{	  
	  const auto px_nextRange = e.first.second/2;
	  const auto idx_nextRange = e.second.second/4;

	  auto add = [&v, px_nextRange, idx_nextRange](auto fromPx, auto fromIdx)
	  {v.emplace_back(std::make_pair(
					 pxR_t(fromPx, px_nextRange),
					 idxR_t(fromIdx, idx_nextRange)));};
	  
	  // top left
	  add(e.first.first, 
	      e.second.first);

	  // top right
	  add(e.first.first+DIM(px_nextRange.x, 0),
	      e.second.first+idx_nextRange);

	  // bottom left
	  add(e.first.first+DIM(0, px_nextRange.y), 
	      e.second.first+2*(idx_nextRange));

	  // bottom right
	  add(e.first.first+px_nextRange,
	      e.second.first+3*idx_nextRange);
	}
      v_levels.insert(v_levels.begin(), v.begin(), v.end());
      std::sort(v_levels.begin(), v_levels.begin()+v.size(),
		// sort with respect to index range
		[](const auto & e0, const auto & e1)
		{return e0.second.first < e1.second.first;});
    }

  hassertm2(helper::ii2n(inDim) == v.size(), helper::ii2n(inDim), v.size());

  std::vector<uint32_t> mp(helper::ii2n(inDim));

  std::for_each(v.begin(), v.end(),
		[&mp, inDim](auto & e)
		{
		  //std::cout << e.first.first << " " << e.first.second << " " << e.second.first << " " << e.second.second << std::endl;
		  assert(e.second.second==1);
		  mp[helper::ii2i(e.first.first, inDim)]=e.second.first;
		}
		);
  
  

  return std::make_tuple(mp, v_levels);
}

template<typename DIM>
auto genMap_grid2qt(DIM inDim)
{
  return std::get<0>(genMap_supergrid2qt(inDim));
}

template<typename QT>
std::pair<uint32_t, uint32_t> getRange_qt2leaves(const uint32_t nodeId, const QT& qt)
{
  const auto level=qt.getLevel(nodeId);
  const auto levelOffset=qt.getLevelOffset(level);

  const auto levelId = nodeId-levelOffset;

  const auto nElems=helper::ipow4(level);
  return std::make_pair(levelId*nElems, nElems);
}


template<typename QT>
auto genMap_qt2leaves(const QT& qt)
{
  std::vector<std::vector<uint32_t>>
    nodes2leaves(qt.nElems());

  for(const auto & leafId :
	helper::range_n(//helper::ii2n(dim_super)
			qt.nElemsLevel_begin
			))
    {
      auto qt_it = qt(leafId);
      do
	{
	  assert(qt_it.idx<nodes2leaves.size());
	  nodes2leaves[qt_it.idx].push_back(leafId);
	  // std::cout << "add leaf " << leafId << " to node "
	  // 	    << qt_it.idx << std::endl;
	}
      while(qt_it());	
    }
  return nodes2leaves;
}

template<typename DIM>
auto gridPos2QTLeaves(const DIM& gridDim, bool swapXY=false)
{
  
  const auto px_idx = std::get<1>(genMap_supergrid2qt(gridDim));;

  std::vector<size_t> o(helper::ii2n(gridDim));
  for(const auto & qtIdx : helper::range_n(o.size()))
    {
      auto e= std::get<0>(px_idx[qtIdx]);
      hassertm(std::get<1>(e).x==1, std::get<1>(e));

      auto pos = std::get<0>(e);

      if(swapXY)
	std::swap(pos.x, pos.y);
      
      o[helper::ii2i(pos, gridDim)] = qtIdx;
      //std::cout << "gridPosId " << std::get<0>(e) << ": qtIdx " << qtIdx  << std::endl;
    }
  return o;
}

template<typename DIM, typename I2>
auto getNeighborsLeavesQT(const DIM dim_super, const std::vector<I2>& off, bool compactify=true)
{
  constexpr bool wrap = false;

  auto inRange = [dim_super](const auto& gridId_o)
  {
    return
      std::min(gridId_o.x, gridId_o.y) >= 0
      && gridId_o.x < dim_super.x
      && gridId_o.y < dim_super.y;
  };
  
  const size_t nNeighbors=off.size();
  std::vector<uint32_t> qtNeighbors(helper::ii2n(dim_super)*nNeighbors, -1);
  {
	
    const auto map_grid2qt = genMap_grid2qt(dim_super);
    for(const auto & gridId : helper::range2_n<int64_t>(dim_super))
      {
	const auto gridIdx = helper::ii2i(gridId, dim_super);
	const auto qtIdx = map_grid2qt[gridIdx];

	size_t cnt=0;

	//std::cout << "qtIdx " << qtIdx  << " at " << gridId << " has these neighbors: ";
	for(const auto& o : off)
	  {
	    auto gridId_o = gridId+o;
	    //std::cout << "[" << gridId_o << "]";

	    auto addNeighbor = [&]()
	    {
	      const auto gridIdx_o = helper::ii2i(gridId_o, dim_super);
	      qtNeighbors[nNeighbors*qtIdx+cnt] = map_grid2qt[gridIdx_o];
	      cnt++;
	    };
	    
	    if(inRange(gridId_o))
	      {
		addNeighbor();		
	      }
	    else if(wrap)
	      {
		if(gridId_o.x >=static_cast<int64_t>(dim_super.x))
		  gridId_o.x-=static_cast<int64_t>(dim_super.x);
		else if(gridId_o.x <static_cast<int64_t>(0))
		  gridId_o.x+=static_cast<int64_t>(dim_super.x);
		
		if(gridId_o.y >=static_cast<int64_t>(dim_super.y))
		  gridId_o.y-=static_cast<int64_t>(dim_super.y);
		else if(gridId_o.y <static_cast<int64_t>(0))
		  gridId_o.y+=static_cast<int64_t>(dim_super.y);
		
		assert(inRange(gridId_o));

		addNeighbor();
	      }
	    else if(!compactify)
	      cnt++;
	  }
	//std::cout << std::endl;
      }
    
  }
  return qtNeighbors;
}

template<typename DIM, typename I2>
auto getNeighborsQT(DIM dim_super, const std::vector<I2>& off, bool compactify=true)
{
  std::vector<std::vector<uint32_t>> o;
  while(dim_super.x >= 1 && dim_super.y >= 1)
    {
      o.push_back(getNeighborsLeavesQT(dim_super, off, compactify));
      dim_super.x/=2;
      dim_super.y/=2;
    }
  return o;
}

template<typename DIM>
auto getParentQT(DIM dim_super)
{
  supertiles_QuadTree<int64_t> qt(dim_super);

  using e_t = uint32_t;
  std::vector<e_t> parentv(qt.nElems(),-1);

  for(const auto & leafId : helper::range_n(qt.nLeaves()))
    {
      auto qt_it=qt(leafId);

      auto childId=leafId;
      while(qt_it() && (parentv[childId]==static_cast<e_t>(-1)))
	{
	  parentv[childId]=qt_it.idx;
	  childId=qt_it.idx;
	}
    }
  return parentv;
}


#endif //__SUPERTILES_QUAD_TREE__
