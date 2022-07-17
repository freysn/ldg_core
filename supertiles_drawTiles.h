#ifndef __SUPERTILES_DRAW_TILES__
#define __SUPERTILES_DRAW_TILES__

#include "supertiles_QuadTree.h"

template<typename T>
V4<uint8_t> x2cairo4(const T& v)
{
  assert(false);
  return V4<uint8_t>(0,0,0,0);
}

template<>
V4<uint8_t> x2cairo4(const double& v)
{
  using CairoDraw_t = helper::CairoDraw<V2<double>>;
  return CairoDraw_t::rgbaf2cairo4<V4<uint8_t>>(V4<double>(v,v,v,1));

  
}

template<typename POS, typename TILE, typename I, typename DIM, typename MAP=std::vector<std::vector<uint32_t>>>
void drawTiles(cairo_t* cr,		
	       TILE* tilesFull,
	       size_t levelOffset,
	       I* tileOffsets,
	       DIM tilePixelDim,
	       uint64_t nTiles,
	       POS* tilesPos,
	       double pixelSize,
	       POS offset,
	       double scale,
	       const MAP& nodes2leaves=MAP())
{
  
  /*
    TILE OFFSETS TO MAP GRID POSITION TO QT IDS
    (POTENTIAL PROBLEM: MAP ONLY WORKS FOR SELECTED LEVEL???)
    USE QT IDS TO MAP DOWN TO LEAVE LEVEL
    CONSIDER THIS
  */

  // for(size_t i=0; i<assignments->size(); i++)
  //   {
  //     std::cout << i;
  //     for(auto it = assignments->cbegin(i);
  // 	  it != assignments->cend(i); it++)
  // 	{
  // 	  const auto supertileId = assignments->getIdx(*it);
  // 	  std::cout << " [" << supertileId << ":" << assignments->getMass(*it) << "]";
  // 	}
  //     std::cout << std::endl;
  //   }

  const bool imgMode=(tilePixelDim.y > 1);

  using CairoDraw_t = helper::CairoDraw<V2<double>>;

  const size_t nPixels = helper::ii2n(tilePixelDim);
  std::vector<V4<uint8_t>> buf;

  if(imgMode)
    buf.resize(nPixels);

  for(const auto i : helper::range_n(nTiles))
    {
      const auto qtIdx = levelOffset+tileOffsets[i];
      const auto pos = offset+scale*tilesPos[i];
      if(imgMode)
	{
	  //const auto tiles=tilesFull+levelOffset*nPixels;
	  
	  for(size_t j=0; j<nPixels; j++)
	    buf[j] =
	      //CairoDraw_t::rgbaf2cairo4<V4<uint8_t>>
	      x2cairo4
	      (tilesFull[qtIdx*nPixels+j]);
	  
	  CairoDraw_t::drawImg(cr, CAIRO_FORMAT_ARGB32,
			       &buf[0],
			       tilePixelDim,
			       pos,
			       scale*pixelSize);
	}
      else
	{
	  // supertiles_DrawSignalCairo<TILE*, DIM> drawSignal(tilesFull, tilePixelDim);
	  // drawSignal.cr = cr;

	  // drawSignal(qtIdx, tilesPos, nodes2leaves);


	  //
	  // draw box
	  //
	  V2<double> chartBoxDim;
	  chartBoxDim.x = chartBoxDim.y = scale*pixelSize*tilePixelDim.x;
	  // signal mode
	  cairo_set_line_width(cr,chartBoxDim.x/100.);
	  cairo_rectangle(cr, pos.x, pos.y, chartBoxDim.x, chartBoxDim.y);	  
	  cairo_set_source_rgba(cr, .0, .0, .0, 1.);
	  cairo_stroke(cr);

	  auto drawSignal = [&](const auto qtIdx)
	  {
	    //
	    // draw signal
	    //
	    cairo_set_source_rgba(cr, .1, .1, .4, 1.);
	    cairo_move_to(cr, pos.x, pos.y+chartBoxDim.y);
	    for(const auto & j : helper::range_n(tilePixelDim.x))
	      cairo_line_to(cr,
			    pos.x+(chartBoxDim.x*j)/tilePixelDim.x,
			    pos.y+((1.-tilesFull[qtIdx*nPixels+j]/*(static_cast<double>(j)/tilePixelDim.x)*/)*chartBoxDim.y)
			    );
	    cairo_stroke(cr);
	  };

	  
	  
	  if(!nodes2leaves.empty())
	    {

	      hassertm3(levelOffset+i < nodes2leaves.size(), levelOffset, i, nodes2leaves.size());
	      const auto leaves = nodes2leaves[qtIdx];
	      std::cout << "there are " << leaves.size() << " leaves\n";
	      for(const auto l : leaves)
		drawSignal(l);
	    }
	  else
	    {
	      drawSignal(qtIdx);
	    }
	}
    }
}

template<typename QT, typename DATA, typename DIM0, typename DIM1, typename N2L, typename PO>
void drawLevels(const QT& qt, const DATA& supertileData, DIM0 dim, DIM1 tileDim, const N2L& nodes2leaves, const PO& po)
{    
  const auto map_qt2grid_leaves =
    helper::invertMap(genMap_grid2qt(dim));
    
  auto qt_it = qt(0);


  size_t cnt=0;

  using CairoDraw_t = helper::CairoDraw<V2<double>>;
    
  do      
    {    

      const auto map_grid2qt = genMap_grid2qt(dim);
      //
      // draw map for debugging
      //
      {
	const double scale = 20.;
	CairoDraw_t cd(helper::cairoBackend_rec);

	auto cr = cd.get();
	cairo_set_font_size(cr, 0.5*scale);
	  
	for(const auto & e : helper::range2_n(dim))
	  {	      
	    cairo_move_to(cr, scale*e.x, scale*e.y);
	    const auto idx_qt = map_grid2qt[helper::ii2i(e, dim)];
	    cairo_show_text(cr, std::to_string(idx_qt).c_str());
	  }

	cd.writePDF(po.outDir+"map_grid2qt_"+std::to_string(cnt)+".pdf");
      }
	
	
      //const auto map_qt2grid = helper::invertMap(map_grid2qt);
	
      CairoDraw_t cd(helper::cairoBackend_rec);

      V2<double> offset(0., 0.);
      std::vector<V2<double>> posv;
      std::vector<size_t> qtIdxv;
      posv.reserve(helper::ii2n(dim));
      qtIdxv.reserve(helper::ii2n(dim));
	

      for(size_t idx_grid : helper::range_n(helper::ii2n(dim)))
	{

	  const auto idx_qt = map_grid2qt[idx_grid];

	  const auto id_grid = helper::i2ii(idx_grid, dim);
	    
	  posv.emplace_back((id_grid.x+.5)/dim.x, (id_grid.y+.5)/dim.y);	    
	  qtIdxv.push_back(idx_qt);	    
	}
	
	
      offset.x=0;
      offset.y=0;
      hassertm3(qt_it.nElemsLevel == helper::ii2n(dim), cnt, qt_it.nElemsLevel, helper::ii2n(dim));	

	
	
      drawTiles(cd.get(), &supertileData[0],
		qt_it.levelOffset,
		&qtIdxv[0], tileDim, posv.size(), &posv[0], 1./(1.1*dim.x*tileDim.x),
		offset, 10.*dim.x, nodes2leaves);
	
      const std::string fname(po.outDir+"supertiles_"+std::to_string(cnt)+".pdf");
      cd.writePDF(fname);

      dim/=2;
      cnt++;
    }
  while(qt_it());
}

#endif //__SUPERTILES_DRAW_TILES__
