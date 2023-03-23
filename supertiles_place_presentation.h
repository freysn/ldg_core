#ifndef __SUPERTILES_PLACE_PRESENTATION__
#define __SUPERTILES_PLACE_PRESENTATION__

#include "supertiles_place_cost.h"
#include "supertiles_QuadTree.h"
#include <set>
#include "helper_CairoDraw.h"
#include "helper_color/cm_plasma.h"
#include "helper_color/cm_map.h"
#include "helper_string.h"
#include "supertiles_place_useCases.h"
#include <map>
#include "supertiles_place_DrawTile.h"
#include "supertiles_place_plan.h"
#include "helper_color.h"
#include "supertiles_place_DrawOpts.h"
#include "supertiles_place_LevelIndicator.h"

namespace supertiles
{
  namespace place
  {

    template<typename DIM1, typename BO>
    auto tileDimDraw(const DIM1& tileDimRef,
		     const uint32_t level,
		     const BO& borderOffsets)
    {
      const auto levelScale=helper::ipow2(level);

      const DIM1 gridPosId(0,0);
      const auto
	borderDelta=DIM1
	(borderOffsets.x[gridPosId.x+levelScale-1]
	 -borderOffsets.x[gridPosId.x],
	 borderOffsets.y[gridPosId.y+levelScale-1]
	 -borderOffsets.y[gridPosId.y]);

      hassertm4(borderDelta.x==borderDelta.y, borderDelta, gridPosId, level, levelScale);
      const auto imgDim =
	V2<double>(borderDelta+tileDimRef*levelScale);

      return imgDim;
    }

    auto level0BorderDefault(uint32_t tileDimX)
    {
     return std::max(tileDimX/20,
	       (uint32_t)1);
    }

    template<typename T0, typename T1>
    auto borderWidthFromLevel(const T0 level, const T1 level0Border)
    {
      //return std::pow(level+1, 1.6)*level0Border;
      return (level+1)*level0Border;
    };

    template<typename DIM, typename B>
    V2<std::vector<uint32_t>> determineBorderOffsets(DIM gridDim, const B level0Border, bool uniformBorder=false)
    {
      V2<std::vector<uint32_t>> borderOffsets;
      auto borderLevel = [](uint32_t tileIdInDim)
      {
	const uint32_t i = tileIdInDim+1;
	// determines lowest set bit
	const uint32_t exp= (i&-i);

	return helper::ilog2(exp);
      };

      if(false)
	{
	  for(const auto i : helper::range_n(16))
	    std::cout << i << " & ";
	  std::cout<<std::endl;
	  for(const auto i : helper::range_n(16))
	    std::cout << borderLevel(i) << " & ";
	  std::cout<<std::endl;
	  exit(0);
	}



      auto borderWidth = [borderLevel, level0Border, uniformBorder](uint32_t tileIdInDim)
      {
	const auto level = borderLevel(tileIdInDim);
	return borderWidthFromLevel(uniformBorder ? 0 : level, level0Border);
      };


      for(const uint32_t dimId : {0,1})
	{

	  uint32_t sum=0;
	  for(const uint32_t tileIdInDim :
		helper::range_n<std::vector<uint32_t>>(gridDim[dimId]))
	    {
	      if(tileIdInDim > 0)
		sum+=borderWidth(tileIdInDim-1);
	      borderOffsets[dimId].push_back(sum);
	    }
	}

      // for(const uint32_t dimId : {0,1})
      //   {
      //     for(const auto e : borderOffsets[dimId])
      //       {
      // 	std::cout << e << " ";
      //       }
      //     std::cout << std::endl;
      //   }
      return borderOffsets;
    }

    template<typename T, typename DIM, typename QT, typename DIM2, typename OV>
    auto determine_leaf2gridPos(const T level, const DIM tileDimRef, const QT qt, const DIM2 gridDim, const OV& borderOffsets, bool swapXY=false)
    {

      

      const auto levelScale=helper::ipow2(level);
      const auto levelDim = gridDim/levelScale;

      const auto leaf2gridPos =
	helper::invertMap(gridPos2QTLeaves(levelDim, swapXY));
      return leaf2gridPos;
    }

      template<typename T, typename DIM, typename QT, typename DIM2, typename OV, typename L2P>
      auto nodeTilePos(const T nodeId, const DIM tileDimRef, const QT qt, const DIM2 gridDim, const OV& borderOffsets, const std::vector<L2P>& leaf2gridPos, bool swapXY=false)	
      {
	const auto level = qt.getLevel(nodeId);

	const auto levelScale=helper::ipow2(level);
	const auto levelDim = gridDim/levelScale;
	
      const auto levelOffset = qt.getLevelOffset(level);
      assert(nodeId >= levelOffset);
      const auto levelNodeId=nodeId-levelOffset;

      assert(levelNodeId < leaf2gridPos.size());

      const auto levelGridPosIdx = leaf2gridPos[levelNodeId];

      const auto gridPosId = helper::i2ii(levelGridPosIdx, levelDim)*levelScale;


      const auto imgPos =
	V2<double>(borderOffsets.x[gridPosId.x]+gridPosId.x*tileDimRef.x,
		   borderOffsets.y[gridPosId.y]+gridPosId.y*tileDimRef.y);

      return imgPos;
    };

    template<typename T, typename DIM, typename QT, typename DIM2, typename OV>
    auto nodeTilePos(const T nodeId, const DIM tileDimRef, const QT qt, const DIM2 gridDim, const OV& borderOffsets, bool swapXY=false)
    {
      const auto level = qt.getLevel(nodeId);
      const auto leaf2gridPos=determine_leaf2gridPos(level, tileDimRef, qt,  gridDim, borderOffsets, swapXY);
      return nodeTilePos(nodeId, tileDimRef, qt, gridDim, borderOffsets, leaf2gridPos, swapXY);
    }

    template<typename nodeId_t, typename disp_t, typename QT, typename IS_VOID>
    auto disparity2activeNodeIndices(const disp_t disparityThreshold,
				     const std::vector<disp_t>& disparities,
				     const std::vector<bool>& drawNode,
				     const QT& qt,
				     const IS_VOID isVoid,
				     const uint32_t levelThresholdMin=-1,
				     const uint32_t levelThresholdMax=-1)
{
  const auto nElemsQT=qt.nElems();
  std::vector<nodeId_t> nodeIndices;
  {
    std::set<nodeId_t> nodeIndicesSet;
    std::vector<bool> splitMask(nElemsQT, false);

    for(const auto & leafId : helper::range_n(qt.nLeaves()))
      {
	if(isVoid(leafId))
	  continue;

	bool split=false;
	auto qt_it = qt(leafId);
	do
	  {
	    if(splitMask[qt_it.idx])
	      break;

	    const auto level = qt.getLevel(qt_it.idx);
	    bool splitHere =
	      (disparities[qt_it.idx] > disparityThreshold)
	      || (level > levelThresholdMax);


	    splitHere = splitHere &&
	      (levelThresholdMin
	       ==
	       static_cast<decltype(levelThresholdMin)>
	       (-1)
	       ||
	       level > levelThresholdMin);

	    split = split || splitHere;

	    splitMask[qt_it.idx]=split;
	  }
	while(qt_it());
      }

    for(const auto & leafId : helper::range_n(qt.nLeaves()))
      {
	if(isVoid(leafId))
	  continue;

	uint32_t prevId = leafId;
	hassertm(!splitMask[leafId]/* || !drawNode[leafId]*/, leafId);

	auto qt_it = qt(leafId);
	do
	  {
	    if(splitMask[qt_it.idx])
	      break;

	    if(!drawNode[qt_it.idx])
	      {
		hassertm(qt_it.idx!=leafId || isVoid(leafId), leafId);
		break;
	      }

	    prevId=qt_it.idx;
	  }
	while(qt_it());
	nodeIndicesSet.insert(prevId);
      }

    nodeIndices=std::vector<nodeId_t>(nodeIndicesSet.begin(), nodeIndicesSet.end());
  }
  return nodeIndices;
}


    double disparityIndicatorLineWidth(const uint32_t level, const double level0Border)
    {
      const double borderWidthFac=.15;
      return borderWidthFromLevel(level, level0Border)*borderWidthFac;
    }

    // template<typename ID, typename QT,
    //   typename DIM, typename DIM2, typename O>
    // auto node2tileImg(ID nodeId, const QT& qt, const double scale,
    // 		       const DIM2 tileDimRef,
    // 		       const DIM& gridDim,
    // 		       const O& borderOffsets,
    // 		       const double level0Border,
    // 		       bool swapXY=false)
    // {

    //   const auto level=qt.getLevel(nodeId);
    //   const auto imgPos = scale*nodeTilePos(nodeId, tileDimRef, qt, gridDim, borderOffsets, swapXY);


    //   const double borderWidthFac=.15;
    //   const auto lineWidth =
    // 	borderWidthFromLevel(level, level0Border)*scale*borderWidthFac;
    //   //(level+1)*level0Border*scale*borderWidthFac;

    //   return std::make_tuple(imgPos, lineWidth);
    // }

    

    template<typename ID, typename PV,
	     typename DISP,
	     typename QT,
	     typename DIM,
      typename DIM2,
      typename O
      //, typename LID
      >
    auto drawDisparityIndicatorNode
    (
     const ID nodeId,

     const PV& parentv,
     const DISP& disparities,


     const QT& qt, //const double scale,
     const DIM2& tileDimRef,
     const DIM& gridDim,
     const O& borderOffsets,
     const double level0Border,
     //const LID& levelImgDim,
     bool swapXY=false)
    {
	const double scale=1.;

	using F2=V2<double>;
	if(parentv[nodeId]==static_cast<decltype(parentv[0])>(-1))
	  return std::make_tuple(F2(0.,0.), static_cast<double>(0),
				 F2(0.,0.), static_cast<double>(0),
				 static_cast<double>(0));

	const auto imgPos = scale*nodeTilePos(nodeId,
					      tileDimRef, qt, gridDim, borderOffsets, swapXY);

	auto imgDim=
	  //std::get<1>(rsltParent);
	  tileDimDraw(tileDimRef,
		      qt.getLevel(nodeId),
		      borderOffsets)*scale;

	const auto parentNodeId = parentv[nodeId];
	const auto imgPosParent = scale*nodeTilePos(parentNodeId,
					      tileDimRef, qt, gridDim, borderOffsets, swapXY);

	const auto lineWidth = scale*disparityIndicatorLineWidth(qt.getLevel(nodeId), level0Border);

	// const auto rsltParent = node2tileImg(parentv[nodeId],
	// 				     qt, scale,
	// 				     tileDimRef,
	// 				      gridDim,
	// 				      borderOffsets,
	// 				      level0Border,
	// 				      swapXY);

	// auto borderDelta=V2<double>(std::get<1>(rsltParent),
	// 				  std::get<1>(rsltParent));

	auto borderDelta=V2<double>(lineWidth,
				    lineWidth);

	borderDelta*=2.;

	//const auto imgPosParent=std::get<0>(rsltParent);
	const auto imgDimParent=
	  //std::get<1>(rsltParent);
	  tileDimDraw(tileDimRef,
		      qt.getLevel(parentNodeId),
		      borderOffsets)*scale;
	const auto imgCenterParent = imgPosParent+.5*imgDimParent;


	const auto imgPosBorderMin=imgPos-borderDelta;
	const auto imgPosBorderMax=imgPos+imgDim+borderDelta;

	//const auto imgPosMax=imgPos+imgDim;
	const auto imgPosMax=imgPosBorderMax;

	//auto originX=imgPos;
	auto originX=imgPosBorderMin;
	originX.y=imgPosBorderMin.y;

	//auto originY=imgPos;
	auto originY=imgPosBorderMin;
	originY.x=imgPosBorderMin.x;

	imgDim=imgDimParent/2;

	uint8_t code=0;
	if(imgPos.x < imgCenterParent.x)
	  {
	    imgDim.x = -imgDim.x;
	    originX.x= imgPosMax.x;
	    originY.x= imgPosBorderMax.x;

	    code |= 1;
	  }



	if(imgPos.y < imgCenterParent.y)
	  {
	    imgDim.y = -imgDim.y;
	    originX.y = imgPosBorderMax.y;
	    originY.y = imgPosMax.y;

	    code |= 2;
	  }

	if(code==2)
	  code=3;
	else if(code==3)
	  code=2;

	//std::cout << "nodeId: " << nodeId << " code " << static_cast<int>(code) <<std::endl;
	const double angleQuart=M_PI/2.;
	double angleOffset=0.;

	angleOffset += code*angleQuart;





	const auto newImgDim = imgDim*disparities[nodeId];

	const auto toX=originX.x+newImgDim.x-
	  (originX.x-imgCenterParent.x);

	const auto toY=originY.y+newImgDim.y-
	  (originY.y-imgCenterParent.y);

	return std::make_tuple(originX, toX, originY, toY, lineWidth);
      }


    

    template<typename ID, typename PV,
	     typename DISP,
	     typename QT,
	     typename DIM,
	     typename DIM2,
	     typename O
	     //, typename LID
	     >
    auto drawLevelIndicatorNode
    (
     const ID nodeId,

     const PV& parentv,
     const DISP& disparities,


     const QT& qt, //const double scale,
     const DIM2& tileDimRef,
     const DIM& gridDim,
     const O& borderOffsets,
     const double level0Border,
     const double levelIndicatorScale,
     //const LID& levelImgDim,
     bool swapXY=false)
    {
      const double scale=1.;
      auto _drawDisparityIndicatorNode = [&](auto n)
      {
	return
	  supertiles::place::drawDisparityIndicatorNode
	  (n,
	   parentv,
	   disparities,
	   qt,
	   tileDimRef,
	   gridDim,
	   borderOffsets,
	   level0Border,
	   //levelImgDim,
	   swapXY);
      };
                  

      const auto parent = qt.parent(nodeId);
      assert(parent != -1);

      const auto n0 = qt.firstChild(parent);
      const auto a = _drawDisparityIndicatorNode(n0);
      const auto b = _drawDisparityIndicatorNode(n0+3);
	
      const auto center= (std::get<0>(a)+std::get<2>(a)+std::get<0>(b)+std::get<2>(b))/4;
      
      const auto level=qt.getLevel(nodeId);

      

      const auto levelIndicatorWH
	= 4*scale*supertiles::place::borderWidthFromLevel(level, level0Border)*levelIndicatorScale;
      
      V4<double> rect(center.x-levelIndicatorWH,
		      center.y-levelIndicatorWH,
		      2*levelIndicatorWH,
		      2*levelIndicatorWH);


      return LevelIndicatorNode(rect, level, qt.nLevels());
    }


      template<typename ID, typename PV, 
      typename DISP,
      typename QT,
      typename DIM,
      typename DIM2,
	       typename O,
	       typename IMG_OPTS
      //,
      //typename LID
      >
      void drawDisparityIndicatorNodeCairo
      (cairo_t* cr,
       const ID nodeId,
       // F2 imgPos,
       // F2 imgDim,
       //F2 borderDelta,
       //const double lineWidth,
       const PV& parentv,
       //const COL& disparityIndicatorCol,
       const DISP& disparities,


       const QT& qt, const double scale,
       const DIM2& tileDimRef,
       const DIM& gridDim,
       const O& borderOffsets,
       const double level0Border,
       //const double levelIndicatorScale,
       //const LID& levelImgDim,
       //bool swapXY/*=false*/,
       const IMG_OPTS& imgOpts)
      {
	const auto swapXY=imgOpts.swapXY;
	const auto levelIndicatorScale=imgOpts.levelIndicatorScale;
	const auto disparityIndicatorCol=imgOpts.disparityIndicatorCol;
	
	const auto [originX, toX, originY, toY, lineWidth]
	  = drawDisparityIndicatorNode
	  (nodeId,
	   // imgPos,
	   // imgDim,
	   //borderDelta,
	   parentv,
	   disparities,
	   qt,
	   //scale,
	   tileDimRef,
	   gridDim,
	   borderOffsets,
	   level0Border,
	   //levelImgDim,
	   swapXY);

	// there is no parent
	if(originX.x == 0 && originX.y == 0
	   && toX==0 &&
	   originY.x==0 && originY.y==0
	   && toY==0)
	  {
	    std::cout << "nodeId " << nodeId << " ... skipping (no parent?)\n";
	    return;
	  }

	cairo_set_line_width(cr, lineWidth*3.5*scale);

	cairo_set_source_rgba(cr,
			      disparityIndicatorCol.x,
			      disparityIndicatorCol.y,
			      disparityIndicatorCol.z,
			      disparityIndicatorCol.w);

	cairo_move_to(cr, originX.x*scale, originX.y*scale);
	cairo_line_to(cr, toX*scale,
		      originX.y*scale);
	cairo_stroke(cr);

	cairo_move_to(cr, originY.x*scale, originY.y*scale);
	cairo_line_to(cr,
		      originY.x*scale,
		      toY*scale);
	cairo_stroke(cr);

	//
	// THE NEW LEVEL INDICATOR
	//
	{
	  const auto lin
	    = drawLevelIndicatorNode
	    (nodeId,
	     // imgPos,
	     // imgDim,
	     //borderDelta,
	     parentv,
	     disparities,
	     qt,
	     //scale,
	     tileDimRef,
	     gridDim,
	     borderOffsets,
	     level0Border,
	     levelIndicatorScale,
	     //levelImgDim,
	     swapXY);

#if 1
	  lin.drawCairo(cr, scale);
#else

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
	  drawTri(lin.root, lin.left, lin.right, lin.lineWidth(), V4<double>(0., 0., 0., 1.));
	  drawTri(lin.root, lin.hleft, lin.hright, lin.hlineWidth(), V4<double>(0.66, 0.66, 0.66, 1.));
	  for(const auto & e : lin.grd)
	    drawLine(std::get<0>(e), std::get<1>(e), lin.grdLineWidth(), V4<double>(0.88, 0.88, 0.88, 1.));
	  
	  // const V2<double> c(originY.x*scale, originY.y*scale);
	  // const auto level = qt.getLevel(nodeId);
	  // const auto w = scale*borderWidthFromLevel(level, level0Border);

	  // cairo_set_source_rgba(cr, 0., 1., 0., 1.);
	  
	  // std::cout << "DRAW LEVEL INDICATOR " << scale << " " <<  c.x-w << " " << c.y-w << " " << 2*w << " LIN " << lin.root << " " << lin.left << " " << lin.right << std::endl;
	  // cairo_rectangle(cr, c.x-w, c.y-w, 2*w, 2*w);
	  // cairo_fill(cr);
#endif
	}
      }    

    template<typename QT>
    auto nodeDrawMask(const QT& qt, const size_t nTilesAssign)
    {
      std::cout << "generate map checking for regular tiles\n";
      const auto isRegular = Plans::isRegularMap<supertiles_QuadTree<int64_t>>(nTilesAssign, qt.nLeaves());
      // std::cout << "generate map checking for regular tiles\n";
      // const auto nodes2leaves=genMap_qt2leaves(qt);

      const auto nNodes=qt.nElems();
      std::vector<bool> drawNode(nNodes, true);

#ifndef NO_OMP
#pragma omp parallel for
#endif
      for(decltype(qt.nElems()) nodeId=0; nodeId<nNodes; nodeId++)
	{
	  //const auto & leaves=nodes2leaves[nodeId];
	  // drawNode[nodeId] = std::all_of(leaves.begin(), leaves.end(),
	  // 				   [&isRegular](size_t n){return isRegular[n];});
	  const auto leafRange=getRange_qt2leaves(nodeId, qt);
	  //for(const auto leafId : leaves)
	  for(auto leafId=leafRange.first; leafId<leafRange.first+leafRange.second; leafId++)
	    {
	      //assert(leafId==leaves[leafId-leafRange.first]);
	      if(!isRegular[leafId])
		{
		  drawNode[nodeId]=false;
		  break;
		}
	    }
	}

      // for(size_t nodeId=0; nodeId<nNodes; nodeId++)
      //   std::cout << "drawNode " << nodeId << ": " << drawNode[nodeId] << std::endl;

      // size_t cnt=0;
      // for(size_t nodeId=0; nodeId<qt.nLeaves(); nodeId++)
      //   cnt+=drawNode[nodeId];

      // std::cout << "LEAFDRAWCOUNT: " << cnt << std::endl;
      return drawNode;
    }

    // use node displarity for color mapping of rectangles; just draw rectangle around tile for indicating deviation in lower levels, draw rectangle around siblings indicating disparity with respect to parent"
    template<typename A, typename DISPARITIES, typename COSTS, typename DIM0, typename DRAW_TILE, typename IMG_OPTS>
    void presentationAdaptiveGrid(
				  const std::string fnamePDF,
				  const std::string fnamePNG,
				  const DIM0 gridDim,
				  const A& qtLeafAssignment,
				  // const std::vector<V4<double>> rep_supertileData,
				  // const REP_OFFV& rep_levelOffsets,
				  // const DIM1 tileDim,
				  const DRAW_TILE& drawTile,
				  const DISPARITIES& disparities,
				  const COSTS& costs,
				  const double disparityThreshold,
				  const uint32_t levelThresholdMin,
				  const uint32_t levelThresholdMax,
				  const double level0Border,
				  const double scale/*=1.*/,
				  const bool drawDisparityIndicator,
				  const std::vector<bool>& drawNode,
				  const IMG_OPTS& imgOpts,
				  size_t shownNodesIdx=-1
				  )
    {

      std::vector<uint8_t> shownNodes;
      if(shownNodesIdx < imgOpts.shownNodes.size())
	shownNodes=imgOpts.shownNodes[shownNodesIdx];
      
#ifdef __NO_OMP
#error "SHOULD NOT BE DEFINED"
#endif
      //#define __NO_OMP

      const bool annotateTiles=imgOpts.annotateTiles;
      const bool writeIndividualTiles=false;
      using nodeId_t = typename DRAW_TILE::nodeId_t;


      //const V4<double> parentIndicatorCol(.36, .15, .20, 1.);

      const auto parentv = getParentQT(gridDim);



      auto isVoid = [&](const auto leafId)
      {
	return qtLeafAssignment[leafId]==voidTileIdx;
      };

      //
      // determine border widths
      //
      V2<std::vector<uint32_t>> borderOffsets;

      borderOffsets=determineBorderOffsets(gridDim, level0Border, imgOpts.do_uniformBorder);

      supertiles_QuadTree<int64_t> qt(gridDim);
      //
      // determine which nodes to show bottom up
      //

      const auto nElemsQT=qt.nElems();

      //std::vector<nodeId_t> nodeIndices;

      std::vector<nodeId_t> nodeIndices;
      if(shownNodes.empty())
	{
	  nodeIndices
	    =disparity2activeNodeIndices<nodeId_t>
	    (disparityThreshold,
	     disparities,
	     drawNode,
	     qt,
	     isVoid,
	     levelThresholdMin,
	     levelThresholdMax);
	}
      else
	{
	  std::cout << "node indices from active nodes\n";
	  for(size_t i=0; i<shownNodes.size(); i++)
	    if(shownNodes[i])
	      nodeIndices.push_back(i);
	}

      if(nodeIndices.size() > imgOpts.maxNTiles)
	{
	  std::cout << "omitting drawing of " << nodeIndices.size() << " tiles, larger than limit of " << imgOpts.maxNTiles << std::endl;
	  return;
	}

      //
      // draw tiles accordingly
      //
      using CairoDraw_t = helper::CairoDraw<V2<double>>;
      CairoDraw_t cd(helper::cairoBackend_rec);

      using cairo4_t = V4<uint8_t>;

      auto cr=cd.get();

      std::vector<V2<double>> levelImgDim;
      for(const auto & level : helper::range_n(qt.nLevels()))
      {
	const auto tileDimRef=drawTile.getTileDim(0);

	const auto imgDim
	  = tileDimDraw(tileDimRef,
			level,
			borderOffsets);

	levelImgDim.push_back(imgDim*scale);
      }

      //
      // write PDF legend of tile size
      //
      if(fnamePDF!="" && imgOpts.do_drawSizeLegend)
	{
	  CairoDraw_t cd(helper::cairoBackend_rec);


	  cairo_set_source_rgba(cd.get(), 0., 0., 0., 1.);

	  const bool dir=0;

	  const auto maxLen=levelImgDim.back()[dir];

	  const auto tickLen=maxLen/40.;

	  cairo_set_line_width(cd.get(), 0.2*tickLen);

	  cairo_set_font_size(cd.get(), tickLen);

	  cairo_move_to(cd.get(), 0,0);
	  cairo_line_to(cd.get(), 0, maxLen);
	  cairo_stroke(cd.get());

	  for(const auto & level : helper::range_n(levelImgDim.size()))
	    {
	      auto nTiles = helper::ipow4(level);
	      cairo_move_to(cd.get(), 0, levelImgDim[level][dir]);
	      cairo_line_to(cd.get(), tickLen, levelImgDim[level][dir]);
	      cairo_stroke(cd.get());

	      cairo_move_to(cd.get(), tickLen*1.05,
			    levelImgDim[level][dir]);


	      std::string postfix;

	      if(nTiles>=1024*1024)
		{
		  nTiles/=(1024*1024);
		  postfix="m";
		}
	      else if(nTiles>=1024)
		{
		  nTiles/=(1024);
		  postfix="k";
		}

	      cairo_show_text(cd.get(),
			      (std::to_string(nTiles)+postfix).c_str());

	      std::cout << "draw tile size legend, level " << level
			<< "l: " << levelImgDim[level][dir]
			<< ", nTiles: " << nTiles << postfix << std::endl;
	    }

	  cd.writePDF(fnamePDF.substr(0, fnamePDF.size()-4)+"_sizes.pdf");
	}

      if(fnamePDF!="" && imgOpts.do_drawCM)
	{
	  const std::string fname=
	    fnamePDF.substr(0, fnamePDF.size()-4)+"_cm";
	  switch(drawTile.getRepAggregationType())
	    {
	    case repAggregationType_mcmc:	    
	      mcmc_cm(fname, imgOpts.mcmc_cm);
	      break;
	    case repAggregationType_chart:
	      stock_cm(fname, cm_turbo, imgOpts);
	      break;
	    default:
	      std::cout << "Warning: no color map output implemented for this case\n";
	    }	  
	}

      auto node2tileImg_ = [&](auto nodeId, int8_t swapXY=-1)
      {
	if(swapXY==-1)
	  swapXY=imgOpts.swapXY;
	
	const auto imgPos = scale*nodeTilePos(nodeId, drawTile.getTileDim(0),
					      qt, gridDim, borderOffsets, swapXY);

	const auto lineWidth = scale*disparityIndicatorLineWidth(qt.getLevel(nodeId), level0Border);
	// const auto [imgPos, lineWidth] =
	//   node2tileImg(nodeId,
	// 		     qt, scale,
	// 		    drawTile.getTileDim(0),
	// 		     gridDim,
	// 		     borderOffsets,
	// 		     level0Border,
	// 		    //levelImgDim,
	// 		     swapXY);

	auto imgDim=levelImgDim[qt.getLevel(nodeId)];

	if(swapXY)
	  std::swap(imgDim.x, imgDim.y);

	return std::make_tuple(imgPos, imgDim,lineWidth, 1.*V2<double>(lineWidth));
      };


      // const auto lineWidth=level0Border*scale*0.5;
      // cairo_set_line_width(cr, lineWidth);

      auto drawDisparityIndicatorNode_ =
	[&](const auto nodeId
	    //,
	    // auto imgPos,
	    // auto imgDim,
	    //auto borderDelta
	    // ,
	    // const auto lineWidth
	    )
	{
	  return drawDisparityIndicatorNodeCairo
	    (cr,
	     nodeId,
	     // imgPos,
	     // imgDim,
	     //borderDelta,
	     //lineWidth,
	     parentv,
	     //imgOpts.disparityIndicatorCol,
	     disparities,
	     qt, scale,
	     drawTile.getTileDim(0),
	     gridDim,
	     borderOffsets,
	     level0Border,
	     //imgOpts.levelIndicatorScale,
	     //imgOpts.swapXY,
	     imgOpts);
	};


      cairo_set_line_cap  (cr, CAIRO_LINE_CAP_ROUND);





      //
      // prevent drawing of rects only consisting of void tiles
      //
      std::vector<bool> fullyVoid(nElemsQT, true);
      std::vector<double> fullVoidRatios(nElemsQT);
      {
	//const auto nodes2leaves=genMap_qt2leaves(qt);


	for(const auto & nodeId : helper::range_n(nElemsQT))
	  {
	    size_t cnt=0;

	    const auto leafRange=getRange_qt2leaves(nodeId, qt);
	    //for(const auto & leafId : nodes2leaves[nodeId])
	    for(auto leafId=leafRange.first; leafId<leafRange.first+leafRange.second; leafId++)
	      if(!isVoid(leafId))
		{
		  cnt++;
		  fullyVoid[nodeId] = false;
		}
	    fullVoidRatios[nodeId] = static_cast<double>(cnt)/
	      /*nodes2leaves[nodeId].size()*/leafRange.second;
	  }
      }

      std::map<uint32_t, std::bitset<4>> parentsChildren;


      struct TileImg
      {
	// using DIM_TILE=
	//   typename std::remove_const<typename std::remove_reference<decltype(drawTile.getTileDim(0))>::type>::type
	//   tileDim;
	std::vector<cairo4_t> buf;
	typename std::remove_const<typename std::remove_reference<decltype(drawTile.getTileDim(0))>::type>::type
	tileDim;
	//DIM_TILE tileDim;
	V2<double> imgPos;
	V2<double> imgDim;
	double tileScale;
	V2<double> borderDelta;
	double lineWidth;

	void free()
	{
	  // trick to actually free up memory from vector
	  std::vector<cairo4_t>().swap(buf);
	}
	
      };

#ifdef __NO_OMP
      const size_t batchSize=1;
#else
      //const size_t batchSize=256;
      const size_t batchSize=nodeIndices.size();
#endif

      std::vector<TileImg> tileImgs(batchSize);

      const std::vector<V2<int64_t>> neighborOffsets
	({V2<int64_t>(-1, 0), V2<int64_t>(1, 0), V2<int64_t>(0, -1), V2<int64_t>(0, 1)});
      //({V2<int64_t>(0, -1), V2<int64_t>(0, -1),V2<int64_t>(0, -1),V2<int64_t>(0, -1)});

      const size_t nNeighbors=neighborOffsets.size();
      
      const auto qtNeighborsPerLevel = getNeighborsQT(gridDim, neighborOffsets, false);
      std::vector<V4<double>> neighborSim(batchSize, V4<double>(0,0,0,0));

      std::vector<bool> nodeIndicesMask(nElemsQT,false);
      for(const auto & e : nodeIndices)
	nodeIndicesMask[e]=true;

      const bool do_experimental_neighborSim=false;
      double maxDisparity=0.;
      if(do_experimental_neighborSim)
      {
	auto ___e = drawTile.get_feat_tileData(0)[0];
	using FE=decltype(___e);
	std::vector<FE> commonRep(drawTile.get_nElemsFeatTile(),0.);
	
	const auto nonVoidLeaves =drawTile.getNonVoidLeaves(qt.nElems()-1);
	for(const auto & leafId : nonVoidLeaves)
	  for(size_t i=0; i<commonRep.size(); i++)
	    commonRep[i]+=drawTile.get_feat_tileData_fromLeaf(leafId)[i];

	constexpr distFuncType_t distFuncType = distFuncType_norm2;
	//distFuncType_cosine_normalized=2,
	using D=double;	      		      
	
	
	
	for(const auto & leafId : nonVoidLeaves)
	  {
	    DistOp<distFuncType, D> dist;
	    for(size_t i=0; i<commonRep.size(); i++)
	      dist(commonRep[i]/nonVoidLeaves.size(),
		   drawTile.get_feat_tileData_fromLeaf(leafId)[i]);
	    maxDisparity+=dist.get();
	  }
	maxDisparity/=nonVoidLeaves.size();
	std::cout << "maxDisparity " << maxDisparity << std::endl;
      }
      
#ifdef __NO_OMP
      std::cout << "start sequential rendering, parallel rendering disabled\n";
#else
      std::cout << "start parallel rendering in batches\n";
#endif
      for(size_t batchIdx=0; batchIdx<nodeIndices.size(); batchIdx+=batchSize)
	//for(const auto & nodeId : nodeIndices)
	{
	  helper::progressBar(batchIdx/static_cast<double>(nodeIndices.size()));

	  const size_t maxIdx=std::min(nodeIndices.size(), batchIdx+batchSize);

	  //
	  // create tile images
	  //
#ifndef __NO_OMP
#pragma omp parallel for
#endif
	  for(size_t idx=batchIdx; idx<maxIdx; idx++)
	    {
	      const auto nodeId=nodeIndices[idx];

	      auto & tileImg=tileImgs[idx-batchIdx];
	      //
	      // draw tile representation
	      //
	      tileImg.tileDim=drawTile.getTileDim(nodeId);
	      const auto nElemsRepTile=helper::ii2n(tileImg.tileDim);

	      //std::cout << "nodeId: " << nodeId << ": nElemsRepTile: " << nElemsRepTile << " " << tileImg.tileDim << std::endl;

	      //size_t levelOffset;

	      std::vector<V4<double>> bufRGBA;
	      bufRGBA.resize(nElemsRepTile);
	      tileImg.buf.resize(nElemsRepTile);

	      std::tie(tileImg.imgPos,
		       tileImg.imgDim, /*levelOffset,*/ tileImg.lineWidth, tileImg.borderDelta)
		= node2tileImg_(nodeId);

	      auto lineWidth=tileImg.lineWidth;

	      cairo_set_line_width(cr, lineWidth);

	      const auto rep_begin=bufRGBA.begin();

	      drawTile(rep_begin, nodeId, shownNodesIdx);

	      // output tile image to file
	      if(imgOpts.individualTileImgs !="")
		{
		  //const std::string fname(fnamePNG.substr(0, fnamePNG.size()-4)+"_"+std::to_string(idx)+".png");

		  const std::string fname
		    (imgOpts.individualTileImgs+"_"+helper::leadingZeros(idx, 6)+".png");
		  helper::cimgWriteNormRGBA(fname,
					    bufRGBA, tileImg.tileDim);

		  tileImg.free();
		  std::cout << "write individual tile image to " << fname << ", skipping rest of drawing ..." << std::endl;
		  continue;
		}

	      //using R = V4<double>;

	      for(const auto & i : helper::range_n(nElemsRepTile))
		tileImg.buf[i] = CairoDraw_t::rgbaf2cairo4<cairo4_t>(rep_begin[i]);




	      const auto tileScale = tileImg.imgDim/V2<double>(tileImg.tileDim);
	      hassertm2(helper::apprEq(tileScale.x,tileScale.y), tileScale.x, tileScale.y);
	      tileImg.tileScale=tileScale.x;


	      //
	      // determine similarity to neighbors
	      //
	      if(do_experimental_neighborSim)
		{
	      const auto level=qt.getLevel(nodeId);

	      auto off = qt.getLevelOffset(level);
	      const auto qt_nElems=qt.nElems();
	      
	      assert(nodeId>=off);
	      
	      //if(idx==0)
	      for(size_t j=0; j<nNeighbors; j++)
		{
		  //using IDX=decltype(qt.nElemsLevel_begin);
		  auto n=qtNeighborsPerLevel[level][(nodeId-off)*nNeighbors+j];
		  //std::cout << "j: " << j << " n: " << n << " " << neighborOffsets[j]<< "  " << off << std::endl;
		  if(n==-1)
		    continue;
		  n+=off;
		  //std::cout << "neighbor " << nodeId << ": " << n << std::endl;
		  hassertm2(n<qt_nElems, n, off);
		  
		  
		  while((n!=-1) && !nodeIndicesMask[n])
		    {
		      hassertm2(n<qt_nElems,n, qt_nElems);
		      n=qt.parent(n);
		    }

		  //std::cout << "test " << nodeId << ": " << n << " " << qt_nElems << std::endl;

		  
		  
		  if(n!=-1 && nodeId < n)
		    {
		      hassertm2(n>=0,n, qt_nElems);
		      const auto a=nodeId;
		      const auto b=n;
		      assert(a>=0 && a<qt_nElems);
		      assert(b>=0 && b<qt_nElems);
		      const auto nElemsFeatTile=drawTile.get_nElemsFeatTile();
		      //const auto a_begin=drawTile.get_feat_tileData(a);
		      //const auto b_begin=drawTile.get_feat_tileData(b);

		      auto ___e = drawTile.get_feat_tileData(0)[0];
		      using FE=decltype(___e);
		      std::vector<FE> commonRep(nElemsFeatTile,0.);

		      const auto a_nonVoidLeaves=drawTile.getNonVoidLeaves(a);
		      const auto b_nonVoidLeaves=drawTile.getNonVoidLeaves(b);
		      
		      auto addToCommonRep = [&](const auto & nonVoidLeaves)
		      {			
			for(const auto & leafId : nonVoidLeaves)
			    for(size_t i=0; i<nElemsFeatTile; i++)
			      commonRep[i]+=drawTile.get_feat_tileData_fromLeaf(leafId)[i];
		      };
		      
		      const size_t nodeLeafCount=
			a_nonVoidLeaves.size()+b_nonVoidLeaves.size();
		      
		      assert(nodeLeafCount);

		      addToCommonRep(a_nonVoidLeaves);
		      addToCommonRep(b_nonVoidLeaves);
		      
		      auto dispPerLeaf = [&](auto it_begin)
		      {	
			using D=double;

			constexpr distFuncType_t distFuncType = distFuncType_norm2;
			//distFuncType_cosine_normalized=2,
		      		      
			DistOp<distFuncType, D> dist;
			
			auto it_ref=commonRep.begin();

			auto it_end=it_begin+nElemsFeatTile;
		      
			for(auto it = it_begin; it!=it_end; it++, it_ref++)
			  dist(*it, (*it_ref)/nodeLeafCount);

			const auto rslt=dist.get();
			hassertm2(std::isfinite(rslt),n, qt.nElems());
			return rslt;
		      };

		      double d=0.;

		      for(const auto & leafId : a_nonVoidLeaves)
			d+=dispPerLeaf(drawTile.get_feat_tileData_fromLeaf(leafId));

		      
		      for(const auto & leafId : b_nonVoidLeaves)
			d+=dispPerLeaf(drawTile.get_feat_tileData(leafId));

		      assert(!a_nonVoidLeaves.empty());
		      assert(!b_nonVoidLeaves.empty());
		      
		      //std::cout << "active larger neighbor for " << nodeId << ": " << n << " (j=" << j << ")" << std::endl;
		      d/=nodeLeafCount;

		      

		      //std::cout << "combined disparity " << d << " " << maxDisparity << " " << disparities[a] << " " << d/maxDisparity << std::endl;

		      d/=maxDisparity;
		      if(d < disparityThreshold)
			neighborSim[idx][j]=d;//1.-d;
		    }
		}
		}
	    }

	  if(imgOpts.individualTileImgs =="")
	    {
	  for(size_t idx=batchIdx; idx<maxIdx; idx++)
	    {
	      const auto tileImgIdx=idx-batchIdx;
	      //const auto nodeId=nodeIndices[idx];
	      //const auto & tileImg=tileImgs[idx-batchIdx];
	      auto tileImg = tileImgs[tileImgIdx];

	      //const auto fullVoidRatio=fullVoidRatios[nodeId];	      
	      //
	      // original tile data prior to modification
	      //
	      const auto dim= tileImgs[tileImgIdx].imgDim;
	      const auto p=tileImgs[tileImgIdx].imgPos;
		      
	      for(size_t j=0; j<neighborOffsets.size(); j++)
		{
		  auto f=neighborSim[idx][j];
		  if(f>0)
		    {
		      f=std::min(f, 1.);
		      //f=1.;
		      //std::cout << "draw j " << j << std::endl;
		      //cairo_set_source_rgba(cr, 1., 0.2, 0.8, 0.3);
		      cairo_set_source_rgba(cr, 0., 0., 0., 0.3);
		      if(j==0)
			cairo_rectangle(cr, p.x, p.y, -dim.x, f*dim.y);
		      else if(j==1)
			cairo_rectangle(cr, p.x+dim.x, p.y, dim.x, f*dim.y);
		      else if(j==2)
			cairo_rectangle(cr, p.x, p.y,f*dim.x, -dim.y);
		      else if(j==3)
			cairo_rectangle(cr, p.x, p.y+dim.y, f*dim.x, dim.y);
		      cairo_fill(cr);
		    }
		}
	    }

	  for(size_t idx=batchIdx; idx<maxIdx; idx++)
	    {
	      const auto tileImgIdx=idx-batchIdx;
	      const auto nodeId=nodeIndices[idx];
	      //const auto & tileImg=tileImgs[idx-batchIdx];
	      auto tileImg = tileImgs[tileImgIdx];


	      auto draw = [&](cairo_t* cr)
	      {

		const auto fullVoidRatio=fullVoidRatios[nodeId];

		if(fullVoidRatio < 1. && imgOpts.fullVoidRatioAdapt)
		  {
		    const auto origSize=tileImg.tileScale*V2<double>(tileImg.tileDim);

		    cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		    cairo_set_line_width(cr, tileImg.lineWidth);
		    cairo_rectangle(cr, tileImg.imgPos.x, tileImg.imgPos.y, origSize.x, origSize.y);
		    cairo_stroke(cr);

		    const auto fac =
		      std::sqrt(fullVoidRatio);

		    

		    

		    const auto newSize= fac*origSize;
		    //recompute offset as asked for by reviewers

		    if(fullVoidRatio>0.)
		      {
			const auto leafRange=getRange_qt2leaves(nodeId, qt);
			//for(const auto & leafId : nodes2leaves[nodeId])

			V2<double> center(0,0);
			{
			  size_t cnt=0;
			  for(auto leafId=leafRange.first; leafId<leafRange.first+leafRange.second; leafId++)
			    if(!isVoid(leafId))
			      {
				const auto posDim=node2tileImg_(leafId);
				center+=std::get<0>(posDim)+0.5*std::get<1>(posDim);
				cnt++;
			      }
			  assert(cnt>0);
			  center/=cnt;
			}

			const auto p=center-newSize*0.5;
			  
			tileImg.imgPos=minv(tileImg.imgPos+origSize-newSize,
					    maxv(tileImg.imgPos,p));
		      }
		    else
		      {
			
			const auto off=(origSize-newSize)*0.5;
			tileImg.imgPos+=off;
		      }
		    tileImg.tileScale*=fac;
		  }

		if(fullVoidRatio > 0.)
		  {

		    assert(tileImg.tileScale > 0);

		    CairoDraw_t::drawImg(cr, CAIRO_FORMAT_ARGB32,
  					 &tileImg.buf[0],
					 tileImg.tileDim,
					 tileImg.imgPos,
					 tileImg.tileScale);
		  }
		if(annotateTiles)
		  {
		    // const auto tileCol=CairoDraw_t::cairo42rgba(tileImg.buf[0]);
		    // const auto tileLuminance=0.3*tileCol.x+0.59*tileCol.y+0.11*tileCol.z;

		    // double fontLuminance=std::sqrt(tileLuminance);
		    // if(tileLuminance<0.5)
		    //   fontLuminance=1.-fontLuminance;

		    const double fontLuminance=1.;

		    const V2<double> dim(tileImg.tileDim.x*tileImg.tileScale,
					 tileImg.tileDim.y*tileImg.tileScale);

		    const auto border=dim*0.08;
		    std::cout << "draw node id " << nodeId << " disparity: " << disparities[nodeId] << std::endl;



		    auto drawText = [&cr, &fontLuminance](const std::string text, const V2<double>& textPos, const size_t fontSize)
		    {
		      cairo_set_font_size(cr, fontSize);
		      #if 0
		      cairo_set_source_rgba(cr, 1,1,1, 0.5);
		      const auto f=std::sqrt(2.)/2.;
		      const auto d=0.05*fontSize;
		      for(const auto delta : {f*V2<double>(d,d), f*V2<double>(d,-d), f*V2<double>(-d,d), f*V2<double>(-d,-d),
					      V2<double>(0,d), V2<double>(0,-d), V2<double>(-d,0), V2<double>(d,0)})
			{
			  cairo_move_to(cr, textPos.x+delta.x, textPos.y+delta.y);
			  cairo_show_text(cr, text.c_str());
			}
		      cairo_set_source_rgb(cr, 0,0,0);
#else
		      cairo_set_source_rgb(cr, fontLuminance,fontLuminance,fontLuminance);
#endif

		      cairo_move_to(cr, textPos.x, textPos.y);
		      cairo_show_text(cr, text.c_str());
		    };

		    const double fontSize=0.5*tileImg.tileScale*tileImg.tileDim.x;
		    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
		    {

		      V2<double> textPos(tileImg.imgPos.x+border.x, tileImg.imgPos.y+fontSize/*+border.y*/);
		      drawText(std::to_string(nodeId), textPos, fontSize);

		    }


		    // cairo_set_font_size(cr, 0.5*fontSize);
		    // cairo_move_to(cr, tileImg.imgPos.x+border.x, tileImg.imgPos.y+dim.y-border.y);

		    std::string s;

		    if(nodeId>=64)
		      {
			s="Δ=";
			std::stringstream ss;
			ss << std::fixed << std::setprecision(2) << disparities[nodeId];
			s+=ss.str();
		      }
		    else if(nodeId<costs.size())
		      {
			assert(nodeId<costs.size());
			std::stringstream ss;
			ss << std::fixed << std::setprecision(2) << costs[nodeId];
			s="γ="+ss.str();
		      }
		    //cairo_show_text(cr, s.c_str());
		    std::cout << "draw string |" << s << "|" << std::endl;
		    drawText(s, V2<double>(tileImg.imgPos.x+border.x, tileImg.imgPos.y+dim.y-border.y), 0.5*fontSize);
		  }
	      };

	      draw(cr);

	      if(writeIndividualTiles)
		{
		  CairoDraw_t cd(helper::cairoBackend_rec);
		  draw(cd.get());

		  //cd.writePDF(fnamePDF.substr(0, fnamePDF.size()-4)+"_"+std::to_string(nodeId)+".pdf");

		  const std::string fname=fnamePDF.substr(0, fnamePDF.size()-4)+"_"+helper::leadingZeros(nodeId, 6)+".pdf";

		  std::cout << "write individual tile image " << fname << std::endl;
		  cd.writePDF(fname);
		}

	      //
	      // draw indicator rectangles: lower level
	      //
	      if(drawDisparityIndicator && imgOpts.disparityIndicatorCol.w > 0.)
		{
		  //auto tileImg = tileImgs[tileImgIdx];
		  drawDisparityIndicatorNode_
		    (nodeId/*, tileImg.imgPos, tileImg.imgDim,*/ /*tileImg.borderDelta, tileImg.lineWidth*/);
		}

	      if(!fullyVoid[nodeId])
		{
		  const auto nodeIdParent = parentv[nodeId];
		  if(nodeIdParent != -1)
		    {
		      auto it = parentsChildren.insert(std::make_pair(nodeIdParent, std::bitset<4>())).first;
		      const auto levelOffset = qt.getLevelOffset(qt.getLevel(nodeId));
		      const auto siblingId = (nodeId-levelOffset)%4;
		      assert(!it->second[siblingId]);
		      it->second.set(siblingId);
		    }
		}
	    }	    
	}
	}

      

#if 0
      //
      // draw parent tile indicator
      //
      for(const auto & parentChildren : parentsChildren)
	{

	  // higher level
	  const auto nodeIdParent = parentChildren.first;
	  const auto childrenNonVoid=parentChildren.second;


	  V2<double> imgPos, imgDim;
	  double lineWidth;

	  std::tie(imgPos, imgDim, std::ignore, std::ignore) = node2tileImg_(nodeIdParent);

	  const auto imgCenter = imgPos+.5*imgDim;

	  imgPos+=
	    0.5*(1-disparities[nodeIdParent])*imgDim;
	  imgDim *= disparities[nodeIdParent];

	  if(drawDisparityIndicator)
	    {
	      cairo_set_line_width(cr, lineWidth);

	      cairo_set_source_rgba(cr,
				    parentIndicatorCol.x,
				    parentIndicatorCol.y,
				    parentIndicatorCol.z,
				    parentIndicatorCol.w);


	      size_t id_00=0;
	      size_t id_10=1;
	      size_t id_01=2;
	      size_t id_11=3;


	      if(swapXY)
		std::swap(id_10, id_01);

	      const bool up=(childrenNonVoid[id_00] || childrenNonVoid[id_10]);
	      const bool down=(childrenNonVoid[id_01] || childrenNonVoid[id_11]);
	      const bool left=(childrenNonVoid[id_00] || childrenNonVoid[id_01]);
	      const bool right=(childrenNonVoid[id_10] || childrenNonVoid[id_11]);

	      if(up)
		{
		  cairo_move_to(cr, /*imgPos.x+0.5*imgDim.x*/imgCenter.x, imgCenter.y);
		  cairo_line_to(cr, /*imgPos.x+0.5*imgDim.x*/imgCenter.x, imgCenter.y-.5*imgDim.y);
		  cairo_stroke(cr);
		}

	      if(down)
		{
		  cairo_move_to(cr, /*imgPos.x+0.5*imgDim.x*/imgCenter.x, imgCenter.y);
		  cairo_line_to(cr, /*imgPos.x+0.5*imgDim.x*/imgCenter.x, imgCenter.y+.5*imgDim.y);
		  cairo_stroke(cr);
		}

	      if(right)
		{
		  cairo_move_to(cr, imgCenter.x, /*imgPos.y+.5*imgDim.y*/imgCenter.y);
		  cairo_line_to(cr, imgCenter.x+.5*imgDim.x, /*imgPos.y+0.5*imgDim.y*/imgCenter.y);
		  cairo_stroke(cr);
		}

	      if(left)
		{
		  cairo_move_to(cr, imgCenter.x, /*imgPos.y+.5*imgDim.y*/imgCenter.y);
		  cairo_line_to(cr, imgCenter.x-.5*imgDim.x, /*imgPos.y+0.5*imgDim.y*/imgCenter.y);
		  cairo_stroke(cr);
		}
	    }
	  //drawRect(nodeIdParent, imgPos-borderDelta, imgDim+2*borderDelta);
	  //}
	}
#endif


      if(imgOpts.individualTileImgs =="")
	{
      // includeInactive=2 : include everything
      auto shownNodesArea_ = [&](const auto& nv, uint32_t
				 includeInactive=false)
      {	
	V2<double> imgPosMin(std::numeric_limits<double>::max(),
			     std::numeric_limits<double>::max());

	V2<double> imgPosMax(std::numeric_limits<double>::lowest(),
			     std::numeric_limits<double>::lowest());

	for(const auto & nodeId : helper::range_n(nv.size()))
	  if(includeInactive==2
	     || (nv[nodeId]==1 || (includeInactive && nv[nodeId]==2)))
	    {
	      V2<double> imgPos, imgDim;
	      std::tie(imgPos,
		       imgDim,
		       std::ignore,
		 std::ignore)
		= node2tileImg_(nodeId, false);

	      // imgPosMin=helper::min(imgPosMin, imgPos);
	      // imgPosMax=helper::max(imgPosMax, imgPos+imgDim);
	      imgPosMin=minv(imgPosMin, imgPos);
	      imgPosMax=maxv(imgPosMax, imgPos+imgDim);
	    }

	if(imgOpts.swapXY)
	  {
	    std::swap(imgPosMin.x, imgPosMin.y);
	    std::swap(imgPosMax.x, imgPosMax.y);
	  }
	return std::make_pair(imgPosMin, imgPosMax-imgPosMin);
      };

      auto shownNodesBorderWidth_=[level0Border, scale, imgOpts]()
      {return level0Border*scale*16.*imgOpts.borderLineWidthScale;};
	

      auto shownNodesBorder_=[&](const auto& nv,
				 const auto& col,
				 bool halfWidth=false)
      {
	const auto [pos, dim]
	  = shownNodesArea_(nv);


	std::cout << nv.size()
		  << " shown nodes, draw rectangle: "
		  << pos << " " << dim
		  << ", col: " << col
		  << std::endl;
	cairo_set_source_rgba(cr,
			      col.x,
			      col.y,
			      col.z,
			      col.w);

	auto lw=shownNodesBorderWidth_();
	if(halfWidth)
	  lw/=2.;
	
	cairo_set_line_width(cr, lw);
	cairo_rectangle(cr,
			pos.x-lw,
			pos.y-lw,
			dim.x+2*lw,
			dim.y+2*lw);
	cairo_stroke(cr);
	return std::make_tuple(pos,dim);
      };



      V2<double> pos(0,0);
      V2<double> dim(0,0);      
      V2<double> border(0,0);
            
      
      if(imgOpts.shownNodes_borderCol.w>0)
	{
	  const auto [p, d]
	    = shownNodesBorder_(shownNodes,
				imgOpts.shownNodes_borderCol);

	  if(imgOpts.zoomActiveNodes)
	    {
	      pos=p;
	      dim=d;
	    }

	  // full line width (as in shownNodesBorder_) and half of if for the middle
	  border=V2<double>(1.5*shownNodesBorderWidth_(),
			    1.5*shownNodesBorderWidth_());
	}
      else if(imgOpts.zoomActiveNodes)
	std::tie(pos,dim) = shownNodesArea_(shownNodes);

      if(imgOpts.zoomActiveNodes)
	{
	  dim.x=std::max(dim.x, dim.y);
	  dim.y=dim.x;
	  
	  V2<double> dimDelta = imgOpts.shownNodes_ghostArea*dim;
	  
	  const auto [p,d] = shownNodesArea_(shownNodes, true);
	  const auto g_from=minv(p+d, maxv(pos-0.5*dimDelta, V2<double>(0., 0.)));
	  const auto g_to=minv(p+d, maxv(pos+dim+0.5*dimDelta,
					 V2<double>(0., 0.)));

	  pos=g_from;
	  dim=g_to-g_from;
	}
      

      {
	const auto & snv=imgOpts.shownNodesv_ann;
	const auto & colv=imgOpts.shownNodesv_ann_borderCol;
	hassertm2(snv.size()==colv.size(), snv.size(), colv.size());

	auto sn_it=snv.begin();
	auto col_it=colv.begin();
	
	for(; sn_it != snv.end() && col_it != colv.end();
	    sn_it++, col_it++)
	{
	  shownNodesBorder_(*sn_it,
			    *col_it,
			    true);
	}
      }

      //
      // add special overlays for agents
      //
      {
	const auto & sn=shownNodes;

	for(size_t nodeId=0; nodeId < sn.size(); nodeId++)
	  {
	    const auto mode=sn[nodeId];
	    if(mode==0)
	      continue;

	    using mark_t=std::uint8_t;
	    enum class MarkMode : mark_t { select = 2 , lookAt=3};
	    

	    V2<double> pos, dim;
	    std::tie(pos,
		     dim,
		     std::ignore,
		     std::ignore)
	      = node2tileImg_(nodeId, false);

	    if(imgOpts.swapXY)
	      {
		std::swap(pos.x, pos.y);
		std::swap(dim.x, dim.y);
	      }

	    switch(mode)
	      {
	      case static_cast<mark_t>(MarkMode::select):
		{
		  std::cout << "draw rect "<< pos << " " << dim << "\n";
		  const V4<double> col(1., 0., 1., 1.);
		  const double lw=dim.x/10.;
		  
		  cairo_set_source_rgba(cr,
					col.x,
					col.y,
					col.z,
					col.w);
		  
		  cairo_set_line_width(cr, lw);
		  cairo_rectangle(cr,
				  pos.x-.5*lw,
				  pos.y-.5*lw,
				  dim.x+lw,
				  dim.y+lw);
		  cairo_stroke(cr);
		}
		break;
	      case static_cast<mark_t>(MarkMode::lookAt):
		{
		  std::cout << "draw circle"<< pos << " " << dim << "\n";
		  const V4<double> col(1., 0., 0., 1.);
		  const double lw=dim.x/10.;
		  
		  cairo_set_source_rgba(cr,
					col.x,
					col.y,
					col.z,
					col.w);
		  
		  cairo_set_line_width(cr, lw);

		  cairo_arc (cr,
			     pos.x+dim.x/2,
			     pos.y+dim.y/2,
			     dim.x/4., 0, 2*M_PI);
		  
		  cairo_stroke(cr);
		}
		break;
	      default:
		break;
	      }
	    
	  }
	
      }
      
      if(fnamePNG != "")
	cd.writePNG(fnamePNG, imgOpts.background,
		    border, pos, dim);

      if(fnamePDF != "")
	cd.writePDF(fnamePDF, imgOpts.background,
		    border, pos, dim);

	}

      //       import cairo

      // img = cairo.ImageSurface.create_from_png("in.png")
      // width = img.get_width()
      // height = img.get_height()

      // imgpat = cairo.SurfacePattern(img)

      // scaler = cairo.Matrix()
      // #1 = 100%; 2 = 50%;0.5 = 200%
      // scaler.scale(2,2) #50% downscale in this case
      // imgpat.set_matrix(scaler)

      // #set resampling filter
      // imgpat.set_filter(cairo.FILTER_BEST)

      // canvas = cairo.ImageSurface(cairo.FORMAT_ARGB32,320,240)
      // ctx = cairo.Context(canvas)

      // ctx.set_source(imgpat)
      // ctx.paint()

      // canvas.write_to_png("out.png")
#ifdef __NO_OMP
#undef __NO_OMP
#endif
    }


    template<distFuncType_t distFuncType,
	     typename DIM_GRID, typename AS, typename FTD,
	     typename RTD,
      typename DIM_TILE, typename DSP,
      typename COST>
    void drawAdaptiveGrids
      (const std::string outDir,
       const DIM_GRID gridDim,
       const AS & qtLeafAssignment,
       const FTD & feat_tileData,
       const size_t nElemsFeatTile,
       const RTD & tileData,
       const DIM_TILE tileDim_in,
       const DSP & nodeDisparities,
       const COST & leafCosts,
       repAggregationType_t repAggregationType,
       double scale,
       const bool drawDisparityIndicator,
       const std::string imgOutOpts,
       const size_t nTilesAssign)
    {

      const supertiles_QuadTree<int64_t> qt(gridDim);

      std::cout << "generate adaptive grid, write results to " << outDir << std::endl;
      
      DrawOpts drawOpts(imgOutOpts);


      const std::vector<bool> drawNode =
	nodeDrawMask(qt, nTilesAssign);

      std::cout << "initialize DrawTile" << std::endl;

      DrawTile<AS, FTD, RTD, DIM_GRID, DIM_TILE> drawTile
	(gridDim, qtLeafAssignment,
	 feat_tileData, nElemsFeatTile,
	 tileData,
	 tileDim_in,
	 repAggregationType,
	 distFuncType,
	 nodeDisparities,
	 drawOpts);

      if(drawOpts.level0Border < 0.)
	drawOpts.level0Border=drawOpts.level0BorderScale*level0BorderDefault((uint32_t)drawTile.getTileDim(0).x);

      if(!drawOpts.shownNodes.empty())
	{
	  for(size_t i=0; i<drawOpts.shownNodes.size(); i++)
	    {
	      std::string fnameBase(outDir/*+"shownNodes"*/);
	      if(drawOpts.shownNodes.size()>1)
		fnameBase+=helper::leadingZeros(i,6);
		
	      // if(swapXY)
	      //   fnameBase+="_swapXY";
	      presentationAdaptiveGrid(drawOpts.do_pdf ? fnameBase+".pdf" : "",
				       drawOpts.do_png ? fnameBase+".png" : "",
				       gridDim,
				       //V2<int64_t>(16,16),
				       qtLeafAssignment,
				       drawTile,
				       nodeDisparities,
				       leafCosts,
				       1.1,
				       -1,
				       std::numeric_limits<uint32_t>::max(),
				       drawOpts.level0Border,
				       scale,
				       drawDisparityIndicator,
				       drawNode,
				       drawOpts,
				       i
				       );

	    }
	  // should not be needed anymore with latest change
	  //drawOpts.shownNodes.clear();
	}


      if(!drawOpts.levels.empty())
      {
	// if(drawOpts.levels.empty())
	//   drawOpts.levels=helper::range_n(qt.nLevels());
	
	// std::reverse(std::begin(drawOpts.levels), std::end(drawOpts.levels));


	for(const auto level : drawOpts.levels)
	  {
	    std::cout << "render level " << level << std::endl;
	    std::string fnameBase(outDir+"adaptiveGridLevels_"+helper::leadingZeros(level, 5));
	    // if(drawOpts.swapXY)
	    //   fnameBase+="_swapXY";

	    presentationAdaptiveGrid(drawOpts.do_pdf ? fnameBase+".pdf" : "",
				     drawOpts.do_png ? fnameBase+".png" : "",
				     gridDim,
				     //V2<int64_t>(16,16),
				     qtLeafAssignment,
				     drawTile,
				     nodeDisparities,
				     leafCosts,
				     1.1,
				     -1,
				     level,
				     drawOpts.level0Border,
				     scale,
				     drawDisparityIndicator,
				     drawNode,
				     drawOpts);
	  }
      }

      if(!drawOpts.disparities.empty())
      {
	// size_t maxStep=20;
	// if(drawOpts.disparityMaxStep!=static_cast<decltype(drawOpts.disparityMaxStep)>(-1))
	//   maxStep=drawOpts.disparityMaxStep;

	// auto steps=helper::range_n(maxStep+1);
	// std::reverse(std::begin(steps), std::end(steps));
	for(const auto disparityThreshold : drawOpts.disparities)
	  {
	    //const double disparityThreshold = disparityThresholdI/static_cast<double>(maxStep);

	    std::string fnameBase(outDir+"adaptiveGrid_"
				  +helper::leadingZeros(int(10000.*disparityThreshold), 5));
	    // if(drawOpts.swapXY)
	    //   fnameBase+="_swapXY";
	    presentationAdaptiveGrid(drawOpts.do_pdf ? fnameBase+".pdf" : "",
				     drawOpts.do_png ? fnameBase+".png" : "",
				     gridDim,
				     //V2<int64_t>(16,16),
				     qtLeafAssignment,
				     drawTile,
				     nodeDisparities,
				     leafCosts,
				     disparityThreshold,
				     drawOpts.levelThresholdMin,
				     std::numeric_limits<uint32_t>::max(),
				     drawOpts.level0Border,
				     scale,
				     drawDisparityIndicator,
				     drawNode,
				     drawOpts);
	  }
      }

    }

  }
}
#endif //__SUPERTILES_PLACE_PRESENTATION__
