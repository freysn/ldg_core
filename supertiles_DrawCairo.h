#ifndef __SUPERTILES_DRAW_CAIRO__
#define __SUPERTILES_DRAW_CAIRO__

#include "helper/helper_CairoDraw.h"
#include "supertiles_DrawBase.h"

// template<typename A, typename POS, typename TILE, typename DIM>
// void tiles_draw(cairo_t* cr,
// 		A* assignments,
// 		TILE* tiles,
// 		POS* tilesPos,
// 		uint64_t nTiles,
// 		DIM tileDim,
// 		TILE* supertiles,
// 		POS* supertilesPos,
// 		uint64_t nSupertiles
// 		)
// {

//   // for(size_t i=0; i<assignments->size(); i++)
//   //   {
//   //     std::cout << i;
//   //     for(auto it = assignments->cbegin(i);
//   // 	  it != assignments->cend(i); it++)
//   // 	{
//   // 	  const auto supertileId = assignments->getIdx(*it);
//   // 	  std::cout << " [" << supertileId << ":" << assignments->getMass(*it) << "]";
//   // 	}
//   //     std::cout << std::endl;
//   //   }

  
  
// }


template<typename DATA, typename DIM>
struct supertiles_DrawSignalCairo : supertiles_DrawBase<DATA,DIM>
{

  template<typename F2, typename F>
  supertiles_DrawSignalCairo(const DATA data, const DIM memberDim,
			     const F2 offset, const F scale)
    :
    supertiles_DrawBase<DATA,DIM>(data, memberDim, offset, scale)
  {
  }

  
  virtual void operator()(uint32_t qtIdx, V2<size_t> tilePos, V2<size_t> memberGridDim, const uint32_t leavesFrom, const uint32_t nLeaves)
  {
    const auto nPixels = helper::ii2n(this->_memberDim);
    
    static_assert(std::is_integral<decltype(tilePos.x)>::value, "Integral required.");

    assert(cr != 0);

    const auto pos = this->grid2canvas(tilePos);

    // supertiles_DrawSignalCairo<TILE*, DIM> drawSignal(tilesFull, memberDim);
    // drawSignal.cr = cr;

    // drawSignal(qtIdx, tilesPos, nodes2leaves);


    //
    // draw box
    //
    V2<double> chartBoxDim;
    chartBoxDim.x = chartBoxDim.y = this->_scale*memberGridDim.x;
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

      cairo_move_to(cr, pos.x, pos.y+chartBoxDim.y);
      for(const auto & j : helper::range_n(this->_memberDim.x))
	cairo_line_to(cr,
		      pos.x+(chartBoxDim.x*j)/this->_memberDim.x,
		      pos.y+((1.-this->_data[qtIdx*nPixels+j]/*(static_cast<double>(j)/_memberDim.x)*/)*chartBoxDim.y)
		      );
      cairo_stroke(cr);
    };

    if(nLeaves > 0)
      {
      cairo_set_source_rgba(cr, .1, .1, 1., 1.);
      for(const auto l : helper::range_bn(leavesFrom, nLeaves))
	drawSignal(l);
      }

    cairo_set_source_rgba(cr, 1., 1., .1, 1.);
    drawSignal(qtIdx);
  }

  cairo_t* cr=0;
};


#endif // __SUPERTILES_DRAW_CAIRO__
