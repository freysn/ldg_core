#ifndef __SUPERTILES_FILTER__
#define __SUPERTILES_FILTER__

#include <numeric>
#include "helper/helper_SAT.h"

//
// TODO: USE THRUST SAT IMPLEMENTATION FOR 2D CASE
//
namespace supertiles
{
  template<typename DATA, typename DIM>
  struct FilterBox
  {
    auto initBase(const DATA& supertileData, const DIM& tileDim_in)
    {
      this->tileDim = tileDim_in;
      this->sums.resize(supertileData.size());

      const auto nElemsPerTile = helper::ii2n(this->tileDim);
      const auto nTiles =
	supertileData.size()
	/nElemsPerTile;
      assert(nTiles*nElemsPerTile==supertileData.size());
      return nTiles;
    }

    template<typename V, typename N, typename R>
    auto outBase(V v, const N nCoveredElems, R radius) const
    {
      assert(nCoveredElems >= 1 && nCoveredElems <= helper::ii2n(tileDim));
      
      v /= nCoveredElems;

      
      auto ratio = static_cast<double>(nCoveredElems);

      if(tileDim.y==1)
	ratio/=(2*radius+1);
      else
	ratio/=(2*radius+1)*(2*radius+1);
      
      assert(ratio>0. && ratio<=1.);
      
      return std::make_tuple(v, ratio);
    }
    
    DATA sums;
    DIM tileDim;
  };

#if 0
  template<typename DATA, typename DIM>
  struct FilterBox1D : FilterBox<DATA,DIM>
  {    
    void init(const DATA& supertileData, const DIM& tileDim_in)
    {
      const auto nTiles = this->initBase(supertileData, tileDim_in);
      for(size_t tileId=0; tileId<nTiles; tileId++)
	std::partial_sum(supertileData.begin()+this->tileDim.x*tileId,
			 supertileData.begin()+this->tileDim.x*(tileId+1),
			 this->sums.begin()+this->tileDim.x*tileId);
    }

    auto operator()(int64_t tileId, int64_t elemId, int64_t radius) const
    {
      //using T = typename std::remove_const<typename std::remove_reference<decltype(this->sums[0])>::type>::type;
      
      const std::pair<int64_t, int64_t> range
	(std::max(elemId-radius-1, static_cast<int64_t>(-1)),
	 std::min(elemId+radius, static_cast<int64_t>(this->tileDim.x-1))
	 );

      const auto s = this->sums.begin()+tileId*this->tileDim.x;
      auto v=s[range.second];

      if(range.first >= 0)
	v-= s[range.first];

      const auto nCoveredElems = range.second-range.first;
      return this->outBase(v, nCoveredElems, radius);
    }
    
  };

#endif

  template<typename DATA, typename DIM>
  struct FilterBox2D : FilterBox<DATA,DIM>
  {    
    void init(const DATA& supertileData, const DIM& tileDim_in)
    {
      const auto nTiles = this->initBase(supertileData, tileDim_in);
      this->sums = supertileData;
      const auto stride = helper::ii2n(this->tileDim);
      
      for(size_t tileId=0; tileId<nTiles; tileId++)
	helper::genSAT(this->sums.begin()+stride*tileId,
		       this->tileDim);
    }

    auto operator()(int64_t tileId, int64_t elemId, int64_t radius) const
    {
      //using T = typename std::remove_const<typename std::remove_reference<decltype(this->sums[0])>::type>::type;

      const V2<int64_t> p = helper::i2ii(elemId, this->tileDim);
      
      const std::pair<V2<int64_t>, V2<int64_t>> range
	(V2<int64_t>(std::max(p.x-radius-1, static_cast<int64_t>(-1)),
		     std::max(p.y-radius-1, static_cast<int64_t>(-1))),
	 V2<int64_t>(std::min(p.x+radius, static_cast<int64_t>(this->tileDim.x-1)),
		     std::min(p.y+radius, static_cast<int64_t>(this->tileDim.y-1)))
	 );

      const auto stride = helper::ii2n(this->tileDim);

      const auto s = this->sums.begin()+tileId*stride;

      const auto v =
	helper::patchSumFromSAT
	(range.first+V2<int64_t>(1,1), range.second+V2<int64_t>(1,1), s, V2<int64_t>(this->tileDim));

      const auto nCoveredElems = helper::ii2n(range.second-range.first);

      return this->outBase(v, nCoveredElems, radius);      
    }
    
  };

  template<typename DATA, typename DIM>
  struct FilterGaussianElem
  {
    auto init(DATA& supertileData, const DIM& tileDim_in)
    {
      tileDim = tileDim_in;
      data = &supertileData;
      rngv.resize(supertileData.size());

      std::mt19937 gen{1337};
      //std::normal_distribution<> dis{v,stddev};
      std::uniform_real_distribution<double> dis(-1., 1.);
      for(auto & e : rngv)
	e = dis(gen);      
    }

    auto operator()(int64_t tileId, int64_t elemId, double radius) const
    {
      // values near the mean are the most likely
      // standard deviation affects the dispersion of generated values from the mean
      

      const auto stride = helper::ii2n(tileDim);

      const auto idx = tileId*stride+elemId;
      const auto v = *(data->begin()+idx);

      assert(v>=0. && v<= 1.);
      //return std::make_tuple(v, /*ratio*/1.);

      //std::random_device rd{};
      //std::mt19937 gen{rd()};

      
      

      //helper::ChronoTimer timer;
      
#if 0
      do
	{
	  v=dis(gen);
	}
      while(v < 0. || v > 1.);
#else
      const auto out=std::max(0., std::min(v+radius*rngv[idx], 1.));
#endif
      // const auto t = timer.get_ms();
      // std::cout << "FILTER took " << timer.get_ms() << " ms" << std::endl;

      assert(radius > 0 || v==out);
      return std::make_tuple(out, /*ratio*/1.);
    }
    
    DATA* data=0;
    DIM tileDim=DIM(0,0);
    std::vector<double> rngv;
  };


  
}

#endif //__SUPERTILES_FILTER__
