#ifndef __SUPERTILES_PLACE_INIT__
#define __SUPERTILES_PLACE_INIT__

namespace supertiles
{
  namespace place
  {
    template<typename NodeLeafCounts, typename D, typename QT, typename A, typename T>
    auto initSupertiles(const QT& qt, const A& qtLeafAssignment, const T& tiles, size_t nElemsTile)
    {
      {
	const double sizeMB=(qt.nElems()*nElemsTile*sizeof(D))*1.e-6;

	if(sizeMB>10000)
	  std::cout << "supertiles rep: allocate " << sizeMB << "MB\n";
      }
      std::vector<D> supertiles(qt.nElems()*nElemsTile, 0.);

      NodeLeafCounts nodeLeafCounts(qt);
      
      //
      // initialize supertiles
      //
      for(const auto & qtIdx : helper::range_n(qt.nLeaves()))
	{
	  auto qt_it = qt(qtIdx);

	  const auto tileIdx = qtLeafAssignment[qtIdx];

	  if(tileIdx != voidTileIdx)
	    {
	      nodeLeafCounts.add(qt_it.idx);
	      do
		{	      
		  for(const auto & elemIdx : helper::range_n(nElemsTile))
		    {
		      auto & out = supertiles[qt_it.idx*nElemsTile+elemIdx];
		      auto & in = tiles[tileIdx*nElemsTile+elemIdx];
		      out += in;
		      //std::cout << "new value at " << qt_it.idx  << ": " << out << " | updated with " << in << " with level offset " << qt_it.levelOffset << std::endl;
		    }
		}
	      while(qt_it());
	    }
	}
      return std::make_tuple(supertiles, nodeLeafCounts);
    }

  }
}

#endif //__SUPERTILES_PLACE_INIT__
