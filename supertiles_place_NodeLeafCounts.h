#ifndef __NODE_LEAF_COUNTS__
#define __NODE_LEAF_COUNTS__

namespace supertiles
{
  namespace place
  {
    //
    // this class can be used for a full tree
    //
    template<typename QT>
    class NodeLeafCounts_full
    {
    public:

      NodeLeafCounts_full(const QT& qt_in=QT(), uint32_t offset_in=0) :
	qt(qt_in), offset(offset_in)
      {}
      
      // gives number of nodes
      auto operator()(uint32_t nodeIdx) const
      {
	const auto level = qt.getLevel(nodeIdx+offset);

	const auto n=QT::nLeavesPerNode(level);
	// if(offset>0)
	//   std::cout << "level: " << level << " " << offset << " " << n << std::endl;
	return n;
      }

      // gives number of nodes
      void add(uint32_t /*leafIdx*/)
      {
	assert(offset==0);
      }

      // gives number of nodes
      void rm(uint32_t /*leafIdx*/)
      {
	assert(offset==0);
      }

      // gives number of nodes
      void addNode(uint32_t /*nodeIdx*/)
      {
	assert(offset==0);
      }

      // gives number of nodes
      void rmNode(uint32_t /*nodeIdx*/)
      {
	assert(offset==0);
      }

#ifndef NDEBUG
      void dbg_sanityCheck()
      {}
#endif

      NodeLeafCounts_full<QT> higherLevelCopy(uint32_t leafLevel) const
      {
	assert(offset==0);
	const auto levelOffset=qt.getLevelOffset(leafLevel);
	assert(qt.getLevel(levelOffset)==leafLevel);
	//std::cout << "LEVEL OFFSET: " << levelOffset << std::endl;
	return NodeLeafCounts_full(qt, levelOffset);
      }

    private:
      //const QT& qt;
      QT qt;
      const uint32_t offset;
    };

    //
    // this class can be used for a tree also featuring empty elements
    //

    
    class NodeLeafCounts_part_copy
    {
    public:
      NodeLeafCounts_part_copy(const int32_t* leafCounts_in) :
	_leafCounts(leafCounts_in)
      {
      }

      // gives number of nodes
      auto operator()(uint32_t nodeIdx) const
      {
	return _leafCounts[nodeIdx];
      }
      
    private:     
      const int32_t* _leafCounts;
    };

    template<typename QT>
    class NodeLeafCounts_part
    {
    public:
      NodeLeafCounts_part(const QT& qt_in) :
	qt(qt_in)
      {
	_leafCounts.resize(qt.nElems(), 0);
      }

      NodeLeafCounts_part()
      {	
      }

      // gives number of nodes
      auto operator()(uint32_t nodeIdx) const
      {
	return _leafCounts[nodeIdx];
      }

      
      void addNode(uint32_t nodeIdx, int32_t delta=1)
      {
	//assert(delta > 0 ||  _leafCounts[nodeIdx] > 0);
	_leafCounts[nodeIdx]+=delta;
	//assert(_leafCounts[nodeIdx] <= QT::nLeavesPerNode(qt.getLevel(nodeIdx)));	
      }

      void rmNode(uint32_t nodeIdx)
      {
	addNode(nodeIdx, -1);
      }
      
      // gives number of nodes
      void add(uint32_t leafIdx, int32_t delta=1)
      {
	auto qt_it = qt(leafIdx);
	do
	  {
	    // assert(delta > 0 ||  _leafCounts[nodeIdx] > 0);
	    // _leafCounts[qt_it.idx]+=delta;
	    // assert(_leafCounts[nodeIdx] <= QT::nLeavesPerNode(qt.getLevel(nodeIdx)));
	    addNode(qt_it.idx, delta);
	  }
	while(qt_it());
      }

      // gives number of nodes
      void rm(uint32_t leafIdx)
      {	
	add(leafIdx, -1);
      }

      NodeLeafCounts_part_copy higherLevelCopy(uint32_t leafLevel) const
      {
	const auto p = &_leafCounts[qt.getLevelOffset(leafLevel)];
	return
	  NodeLeafCounts_part_copy(p);
      }

#ifndef NDEBUG
      void dbg_sanityCheck()
      {}
#endif

      #ifdef NDEBUG
      private:
      #endif
      
      std::vector<int32_t> _leafCounts;
      //const QT& qt;
      QT qt;
    };
  


  
  template<typename QT>
  class DBG_NodeLeafCounts_part_check_copy
  {
  public:
    DBG_NodeLeafCounts_part_check_copy(NodeLeafCounts_part_copy part_copy_in,
				       NodeLeafCounts_full<QT> full_in):
      part_copy(part_copy_in), full(full_in)
    {}
      
    // gives number of nodes
    auto operator()(uint32_t nodeIdx) const
    {
      const auto p = part_copy(nodeIdx);
      const auto f = full(nodeIdx);

#if 1
      hassertm3(p>=0 && p<=f, p, f, nodeIdx);
#else
      hassertm3(p==f, p, f, nodeIdx);
#endif
      
      return p;
    }
    
  private:

    const NodeLeafCounts_part_copy part_copy;
    const NodeLeafCounts_full<QT> full;
  };
  
  template<typename QT>
  class DBG_NodeLeafCounts_part_check
  {
  public:

    DBG_NodeLeafCounts_part_check(const QT& qt_in=QT()) :
      full(qt_in), part(qt_in), qt(qt_in)
    {}

    // gives number of nodes
    auto operator()(uint32_t nodeIdx) const
    {
      const auto p = part(nodeIdx);
      const auto f = full(nodeIdx);
#if 1
      hassertm3(p>=0 && p<=f, p, f, nodeIdx);
#else
      hassertm3(p==f, p, f, nodeIdx);
#endif
      
      return p;
    }

      
    void addNode(uint32_t nodeIdx, int32_t delta=1)
    {
      part.addNode(nodeIdx, delta);
    }

    void rmNode(uint32_t nodeIdx)
    {
      addNode(nodeIdx, -1);
    }
      
    // gives number of nodes
    void add(uint32_t leafIdx, int32_t delta=1)
    {
      part.add(leafIdx, delta);
    }

    // gives number of nodes
    void rm(uint32_t leafIdx)
    {	
      add(leafIdx, -1);
    }

    DBG_NodeLeafCounts_part_check_copy<QT> higherLevelCopy(uint32_t leafLevel) const
    {
      return DBG_NodeLeafCounts_part_check_copy<QT>
	(part.higherLevelCopy(leafLevel),
	 full.higherLevelCopy(leafLevel));
    }

#ifndef NDEBUG
    void dbg_sanityCheck()
    {
      for(const auto leafId : helper::range_n(qt.nElems()))
	(*this)(leafId);
    }
#endif

  private:
    NodeLeafCounts_part<QT> part;
    NodeLeafCounts_full<QT> full;

    //const QT& qt;
    QT qt;
  };
  }
}

#endif //__NODE_LEAF_COUNTS__
