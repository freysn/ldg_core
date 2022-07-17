#ifndef __SUPERTILES_PLACE_LEVEL_INDICATOR_NODE__
#define __SUPERTILES_PLACE_LEVEL_INDICATOR_NODE__
namespace supertiles
{
  namespace place
  {
    struct LevelIndicatorNode
    {
      LevelIndicatorNode() {}

      template<typename R>
      LevelIndicatorNode(const R& rect, const uint32_t level, const uint32_t nLevels)
      {
	init(rect, level, nLevels);
      }
      
      template<typename R>
      void init(const R& rect, const uint32_t level, const uint32_t nLevels)
      {
	if(rect.z==0 || rect.w==0)
	  return;
	
	root=V2<double>(rect.x + (rect.z / 2), rect.y);
	left=V2<double>(rect.x,rect.y+rect.w);
	right=V2<double>(rect.x+rect.z,rect.y+rect.w);
                  

	auto triHeight = [&](auto l)
	{
	  const double s = (nLevels-l-1)/static_cast<double>(nLevels-1);
	  return std::make_tuple(root+s*(left-root),
				 root+s*(right-root));
      
	};

	std::tie(hleft, hright)
	  =
	  triHeight(level);

	for(const auto & l : helper::range_be(1,nLevels-1))
	  if(l!=level)
	    grd.emplace_back(triHeight(l));

	const auto center=V2<double>(rect.x+0.5*rect.z, rect.y+0.5*rect.w);
	const auto actualCenter=(root+left+right)/3;
	const auto delta=center-actualCenter;
      
	trans([delta](auto p){return p+delta;});

	drawLevelLine = (level!=0) && (level!=nLevels-1);
      }
            
      template<typename F>
      void trans(F f)
      {
	root=f(root);
	left=f(left);
	right=f(right);
	hleft=f(hleft);
	hright=f(hright);
	for(auto & e : grd)
	  {
	    std::get<0>(e)=f(std::get<0>(e));
	    std::get<1>(e)=f(std::get<1>(e));
	  }
      }

      double lineWidth() const
      {
	return lineWidthScale*(right.x-left.x)/20.;
      }

      double hlineWidth() const
      {
	//return (hright.x-hleft.x)/20.;
	return lineWidth();
      }

      double grdLineWidth() const
      {
	//return (hright.x-hleft.x)/20.;
	return 0.3*lineWidth();
      }


      void drawCairo(cairo_t* cr, const double scale=1.) const
      {
	if(right==left)
	  return;
	
	auto drawLine=[&](auto l, auto r, auto lw, auto col)
	{	
	  cairo_move_to(cr, l.x*scale, l.y*scale);
	  cairo_line_to(cr, r.x*scale, r.y*scale);

	  cairo_set_source_rgba(cr, col.x, col.y, col.z, col.w);

	  cairo_set_line_width(cr, scale*lw);
	    	    
	  cairo_stroke(cr);
	};
	  
	auto drawTri=[&](auto ro, auto l, auto r, auto lw, auto col)
	{	
	  cairo_move_to(cr, ro.x*scale, ro.y*scale);
	  cairo_line_to(cr, l.x*scale, l.y*scale);
	  cairo_line_to(cr, r.x*scale, r.y*scale);
	  //cairo_move_to(cr, ro.x, ro.y);
	  cairo_close_path(cr);

	  cairo_set_source_rgba(cr, col.x, col.y, col.z, col.w);
	  cairo_fill_preserve(cr);

	  cairo_set_source_rgba(cr, 1., 1., 1., 1.);

	  cairo_set_line_width(cr, scale*lw);
	    	    
	  cairo_stroke(cr);
	};

	cairo_set_line_cap  (cr, CAIRO_LINE_CAP_ROUND);
	// https://www.w3.org/TR/SVG11/types.html#ColorKeywords
	drawTri(root, left, right, /*lineWidth()*/0., V4<double>(0., 0., 0., 1.));
	drawTri(root, hleft, hright, /*hlineWidth()*/0., V4<double>(0.66, 0.66, 0.66, 1.));
	if(drawLevelLine)
	  drawLine(hleft, hright, 1.5*hlineWidth(), V4<double>(1., 1., 1., 1.));
	for(const auto & e : grd)
	  drawLine(std::get<0>(e), std::get<1>(e), grdLineWidth(), /*V4<double>(1., 0., 0., 1.)*/V4<double>(1., 1., 1., 1.));
      }

      V2<double> root, left, right, hleft, hright;
      std::vector<std::tuple<V2<double>, V2<double>>> grd;
      bool drawLevelLine=true;
      double lineWidthScale=1.;
    };
  }
}
#endif
