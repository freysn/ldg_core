#ifndef __SUPERTILES_DIST_FUNC__
#define __SUPERTILES_DIST_FUNC__

#include <unordered_map>
#include "helper/helper_hash.h"
#include "helper/helper_assert.h"
#include "helper/helper_util.h"
#include "supertiles_configTypes.h"

//#define DBG_SUPERTILES_DF

struct supertiles_Hash
{
  std::size_t operator()(const float3& e) const noexcept
  {
    size_t seed=1337;
    helper::hash_combine(seed, e.x);
    helper::hash_combine(seed, e.y);
    #warning "FOR TESTING FOR NOW, EVEN IN FLOAT3 DO NOT CONSIDER THIRD COMPONENT WHEN HASHING"
    if(false)
      helper::hash_combine(seed, e.z);
    return seed;
  }

  std::size_t operator()(const float2& e) const noexcept
  {
    size_t seed=1337;
    helper::hash_combine(seed, e.x);
    helper::hash_combine(seed, e.y);
    return seed;
  }
};

bool operator==(const float2& lhs, const float2& rhs) {
  return lhs.x==rhs.x && lhs.y==rhs.y;
}
bool operator==(const float3& lhs, const float3& rhs) {
  return lhs.x==rhs.x && lhs.y==rhs.y && lhs.z==rhs.z;
}

// template <class T>
// T supertiles_voidPos()
// {
//   return T();
// }

// template<>
// float3 supertiles_voidPos()
// {
//   return make_float3
//     (std::numeric_limits<float>::max(),
//      std::numeric_limits<float>::max(),
//      std::numeric_limits<float>::max());;
// }


  template<distFuncType_t distFuncType, typename T, typename TILE, typename QT>
    class supertiles_DistFunc
    {
      using A = ::intervol::AssignmentsPlain<::intervol::assignment_e_t>;
    public:
      
      typedef T pos_t;
      template<typename DIM>
      supertiles_DistFunc(TILE* tiles, TILE* supertiles, unsigned int* supertilesCnts, unsigned int nElems_tile, const QT& qt_in, const unsigned int* qtNeighbors_in, const DIM tileGridDim) :
	qt(qt_in),
	qtNeighbors(qtNeighbors_in)
      {
	this->pos_source = 0;
	this->pos_target = 0;
	this->nElems_source = 0;
	this->nElems_target = 0;
	this->tiles = tiles;
	this->supertiles = supertiles;
	this->supertilesCnts = supertilesCnts;

	this->nElems_tile = nElems_tile;

	this->tileGridDim.x = tileGridDim.x;
	this->tileGridDim.y = tileGridDim.y;
      }

#ifdef __NVCC__
      __host__ __device__
#endif
      void init(T* pos_source, T* pos_target, unsigned int nElems_source,
		unsigned int nElems_target, T defaultValue, A* assignments)
      {
	std::cout << "assignments: " << *assignments << std::endl;
	this->pos_source = pos_source;
	this->pos_target = pos_target;
	this->nElems_source = nElems_source;
	this->nElems_target  = nElems_target;
	this->defaultValue = defaultValue;
	//this->assignments = assignments;

#ifdef DBG_SUPERTILES_DF
	dbg_supertileAssignments =
	  new std::vector<std::set<uint32_t>>(nElems_target);
#endif	  
	size_t cnt=0;
	auto insert_pos2idx =
	  [&](T& e)
	  {
	    //std::cout << "insert " << e.x << " " << e.y << std::endl;
	    const auto success = pos2idx.insert(std::make_pair(e,cnt)).second;
	    cnt++;
	    if(!success)
	      {
		throw std::runtime_error("identical positions, collision");
	      }	    
	  };
	std::for_each(pos_source, pos_source+nElems_source,
		      insert_pos2idx);
	std::for_each(pos_target, pos_target+nElems_target,
		      insert_pos2idx);
	
	// initialize supertiles
	for(size_t i=0; i<assignments->size(); i++)
	  {
	    std::cout << i;
	    for(auto it = assignments->cbegin(i);
		it != assignments->cend(i); it++)
	      {
		const auto supertileId = assignments->getIdx(*it);
		std::cout << " [" << supertileId << ":" << assignments->getMass(*it) << "]";

		auto qt_it = qt(target2supertileId(supertileId));

		// only do it on the lowest level for debugging purposes
		updateQTLevel(qt_it.idx, i, true);
		
		while(qt_it())
		  {
		    updateQTLevel(qt_it.idx, i, true);
		  }
	      }
	    std::cout << std::endl;
	  }
      }
  
      //  template<typename pos_t>
      double idx_idx(const unsigned int& idx_source, const unsigned int& idx_target) const
      {
	assert(idx_source < nElems_source);
	assert(idx_target >= nElems_source);    
	const T a = getPosSource(idx_source);//pos_source[idx_source];
	const T b = getPosTarget(idx_target);//pos_target[localTargetIdx(idx_target)];
	//return idiffs_distFunc::distFunc(a,b);
	return distOp(a,b, -1, -1);
      }

#if 1
#ifdef __NVCC__
    __host__ __device__
#endif
      double pos_idx(const T& a, const unsigned int& idx_target)
      {
	//return pos_pos(a, pos_target[localTargetIdx(idx_target)]);
	return distOp(getIdxFromPos(a), getIdxFromPos(getPosTarget(idx_target)), -1, -1);
      }
#endif

#ifdef __NVCC__
      __host__ __device__
#endif
      T getPosSource(const unsigned int& idx_source)
      {
	//if(idx_source >= nElems_source)
	//assert(defaultValue == pos_source[idx_source]);
	
	if(idx_source < nElems_source)
	  return pos_source[idx_source];
	else
	  return defaultValue;
      }

#ifdef __NVCC__
      __host__ __device__
#endif
      T* getPosSource_p()
      {
	return pos_source;
      }

      #ifdef __NVCC__
      __host__ __device__
#endif
      const T* getPosSource_pc() const
      {
	return pos_source;
      }

#ifdef __NVCC__
      __host__ __device__
#endif
      T* getPosTarget_p()
      {
	return pos_target;
      }

      #ifdef __NVCC__
      __host__ __device__
#endif
      const T* getPosTarget_pc() const
      {
	return pos_target;
      }

      #ifdef __NVCC__
      __host__ __device__
#endif
      thrust::tuple<T,bool> getPosTarget_inRange(const unsigned int& idx_target)
      {
	unsigned int idx_target_local = localTargetIdx(idx_target);

      //if(idx_target_local >= nElems_target)
      //assert(defaultValue == pos_target[idx_target_local]);
      
      if(idx_target_local < nElems_target)
	return thrust::make_tuple(pos_target[idx_target_local], true);
      else
	return thrust::make_tuple(defaultValue, false);
      }
      
#ifdef __NVCC__
      __host__ __device__
#endif
      T getPosTarget(const unsigned int& idx_target)
      {
	return thrust::get<0>(getPosTarget_inRange(idx_target));
	/*
	unsigned int idx_target_local = localTargetIdx(idx_target);

      //if(idx_target_local >= nElems_target)
      //assert(defaultValue == pos_target[idx_target_local]);
      
      if(idx_target_local < nElems_target)
	  return pos_target[idx_target_local];       
      else
	return defaultValue;
	*/
      }

      unsigned int nElems_source;
      unsigned int nElems_target;

#ifdef __NVCC__
      __host__ __device__
#endif	
      auto getIdxFromPos(const T& e)
      {
	// CURRENTLY, THINGS ARE ONLY TESTED WITH 2D DATA BUT 3D FLOATS ARE USED
	assert(e.z==0);
	
	const auto it = pos2idx.find(e);
	//std::cout << "pos " << e.x << " " << e.y << std::endl;
	hassertm3(it != pos2idx.end(), e.x, e.y, e.z);
	return it->second;
      }


      //template<typename pos_t>
#ifdef __NVCC__
      __host__ __device__
#endif
      double distOp_pos(const T& tileFrom, const T& superFrom, const T& tileTo, const T& superTo)
      {
	//return pos_pos(a,b);
	return distOp(
		      getIdxFromPos(tileFrom),
		      getIdxFromPos(superFrom),
		      getIdxFromPos(tileTo),
		      getIdxFromPos(superTo));
      }
      
      //template<typename pos_t>
#ifdef __NVCC__
      __host__ __device__
#endif
      double operator()(const T& tileFrom, const T& superFrom, const T& tileTo, const T& superTo)
      {	
	return distOp_pos(tileFrom, superFrom, tileTo, superTo);
      }

      
    private:
#ifdef __NVCC__
    __host__ __device__
#endif
      unsigned int localTargetIdx(const unsigned int& idx_target)
      {
#ifndef USE_CUDA_RUN_DEVICE
	using namespace std;
#endif
      //assert(idx_target >= nElems_source);

#ifndef USE_IDIFFS_HACK_IV
	assert(idx_target >= nElems_source);
	return idx_target-nElems_source;
#else
	#warning "THIS SEEMS TO BE REQUIRED FOR IDIFFS, BUT ACTUALLY DOESNT SEEM LIKE A TOO GOOD IDEA IN GENERAL, FIX THAT!!"
	//printf("localTargetIdx %d, %d, target %d", nElems_source, nElems_target, idx_target);
	assert(idx_target >= max(nElems_source, nElems_target));
	return idx_target-
	//nElems_source
	  max(nElems_source, nElems_target)
	;
#endif
      }

      T defaultValue;

      T* pos_source;
      T* pos_target;

      const QT qt;
      
      std::unordered_map<T, uint64_t, supertiles_Hash> pos2idx;

      /*
	IDEA:
	conceptually remove both tiles from supertiles QT
	on this basis, then compute distances
	(this circumvents a certain type of bias toward the current location)
       */
      double distOp(
		    const uint32_t idx_tileFrom,
		    uint32_t idx_supertileFrom,
		    const uint32_t idx_tileTo,
		    uint32_t idx_supertileTo)
      {
	// original tiles have lower ids
	// supertiles larger ids accordingly
	//assert(st::max(idx_tileFrom, idx_tileTo) < std::min(idx_supertileFrom, idx_supertileTo));

	assert(idx_tileFrom < idx_supertileFrom);		
	idx_supertileFrom = localTargetIdx(idx_supertileFrom);

	assert((idx_tileTo==-1)==(idx_supertileTo==-1));

	const bool individualMode=(idx_tileTo==-1);

	if(!individualMode)
	  idx_supertileTo = localTargetIdx(idx_supertileTo);
	
	auto distTILE = [&](unsigned int idx_tile,
			    unsigned int idx_supertile,
			    unsigned idx_supertile_alt,
			    bool considerLeaveLevel,
			    uint32_t idx_tile_removeBeforeComparison,
			    uint32_t idx_supertile_removeBeforeComparisonIf)
	{	  
	  // TRAVERSE SUPERTILE WITH QUADTREE
	  
	  assert(idx_tile < qt.nLeaves());
	  
	  auto qt_it = qt(idx_supertile);

	  // stop when the level is reached where that
	  // both considered tiles are currently assigned to
	  auto qt_it_stop = qt(idx_supertile_alt);

	  auto qt_it_remove = qt(idx_supertile_removeBeforeComparisonIf);

	  double sum=0.;
	  //traverse first, not interested at leave level
	  //while(qt_it())

	  //size_t dbg_levelCnt=0;
	  auto nextLevel = [&qt_it, &qt_it_stop, &qt_it_remove]()
	  {
	    const bool success0 = qt_it();
	    const bool success1 = qt_it_stop();
	    const bool success2 = qt_it_remove();

	    assert(success0==success1 || qt_it_stop.isVoid());
	    assert(success0==success2 || qt_it_remove.isVoid());

	    const bool carryOn = (qt_it.idx != qt_it_stop.idx);

	    return carryOn && success0;
	  };


	  if(considerLeaveLevel || nextLevel())
	    do
	      {
		unsigned int cnt = 0;
		const bool removeBeforeComparison =
		  (idx_tile_removeBeforeComparison != -1)
		  &&
		  (qt_it_remove.idx == -1 || qt_it_remove.idx == qt_it.idx);

		if(removeBeforeComparison)
		  {
		    cnt = supertilesCnts[qt_it.idx];
		    assert(cnt>0);
		  }

		
		if(cnt != 1)
		  {
		    double d=0.;
		    double mag=0.;
#ifndef NDEBUG
		    double mag_tile=0.;
#endif
		    for(uint32_t i=0; i<nElems_tile; i++)
		      {
		  
			auto s  = supertiles[qt_it.idx*nElems_tile+i];
			assert(sanityCheckValue(s));

			helper::hassert_norm(s);
			
			const auto & t = tiles[idx_tile*nElems_tile+i];

			helper::hassert_norm(t);
			
			if(removeBeforeComparison)
			  {
			    //VIRTUALLY REMOVE TILE FROM SUPERTILE BEFORE COMPARING
			    const auto dbg_s_old = s;
			    s = getUpdatedValue(s, cnt, -1,
						tiles[idx_tile_removeBeforeComparison*nElems_tile+i]);
			    hassertm4(sanityCheckValue(s), s , dbg_s_old, t, cnt);
			  }

			helper::hassert_norm(s);
		    
			switch(distFuncType)
			  {
			  case distFuncType_norm2:
			    {
			      const auto v = s-t;
			      const auto dc = dot(v,v);
			      d +=dc;
			      hassertm3(dc<= 3.1, dc, s, t);
			    }
			    break;
			  case distFuncType_cosine_normalized:
			    d += s*t;
			    mag += s*s;
#ifndef NDEBUG
			    mag_tile+=t*t;
#endif
			    break;
			  default:
			    assert(false);
			  }
		    
			// std::cout << "super: "<< s.x << " " << s.y << " " << s.z << " " << s.w << " | " << untilSame << std::endl;
			// std::cout << "tile: " << t.x << " " << t.y << " " << t.z << " " << t.w << std::endl;
			// std::cout << "v: " << v.x << " " << v.y << " " << v.z << " " << v.w << std::endl;
		    
		    
		      }
		    if(distFuncType==distFuncType_norm2)
		      d/=nElems_tile;
		    else if(distFuncType==distFuncType_cosine_normalized)
		      {			

#ifndef NDEBUG
			hassertm(mag_tile>=1.-1e-5 && mag_tile <= 1.+1e-5, mag_tile);
#endif

			if(removeBeforeComparison)
			  d/=std::sqrt(mag);
			hassertm(d>=0. && d<=1., d);
			d=1.-d;
		      }
		    sum += d;
		  }
	      }
	    while(nextLevel());
	  	  
	  return sum;
	};
	
	auto distOp_thisPos = [&](auto idx_consideredTile)
	{
	  return distTILE(idx_consideredTile, idx_supertileFrom, idx_supertileTo, false, idx_tileFrom, -1);
	};

	double d=0.;
	
	if(individualMode)
	  d = distOp_thisPos(idx_tileFrom);
	else
	  {
	  d = 
	    distOp_thisPos(idx_tileTo)
	    -
	    distOp_thisPos(idx_tileFrom);
	  //std::cout << "d ( " << idx_tileFrom << " " << idx_tileTo << "): " << d << std::endl;
	  }
	
	// consider neighboring (super)tiles
	if(true)
	  {
	    //check that only supertiles are considered that currently considered tile is NOT a part of
	    const auto otherFac = 0.25;
	    for(uint32_t i=0; i<4; i++)
	      {
		const auto neighbor = qtNeighbors[4*idx_supertileFrom+i];
		if(neighbor != -1)
		  {

		    auto distOp_neighbor = [&](auto idx_consideredTile)
		    {
		      return distTILE(idx_consideredTile, neighbor, idx_supertileFrom,
				      true, idx_tileTo, idx_supertileTo);
		    };

		    if(individualMode)
		      {
			d+=distOp_neighbor(idx_tileFrom);
		      }
		    else
		      d += otherFac*
			(
			 distOp_neighbor(idx_tileTo)
			 -distOp_neighbor(idx_tileFrom)
			 );
		  }
	      }
	  }
	return d;
      }

      TILE getUpdatedValue(const TILE& v_super,
			   const unsigned int oldCnt,
			   const int mod,
			   const TILE& v_tile)
      {
	const auto newCnt = oldCnt+mod;
	if(newCnt==0)
	  return 0;
	else
	  return (oldCnt*v_super+mod*v_tile)/newCnt;
      }


      static bool sanityCheckValue(float4 v)
      {
	const float eps_f = 1e-5;
	return
	  v.x >= -eps_f && v.x <= 1.+eps_f &&
	  v.y >= -eps_f && v.y <= 1.+eps_f &&
	  v.z >= -eps_f && v.z <= 1.+eps_f &&
	  ((v.w >= 1.-eps_f && v.z <= 1.+eps_f)
	   ||(v.w >=-eps_f && v.z <= eps_f)
	   );
      }

      static bool sanityCheckValue(double v)
      {
	//const float eps_f = 1e-5;
	const double eps_f = 1e-6;
	return
	  v >= -eps_f && v <= 1.+eps_f;	   
      }

      void updateQTLevel(uint32_t supertileId,
			   uint32_t tileId,
			   bool add)
      {		

	hassertm2(supertileId < qt.nElems(), supertileId, qt.nElems());

#ifdef DBG_SUPERTILES_DF
	assert((*dbg_supertileAssignments)[supertileId].size()==supertilesCnts[supertileId]);

	auto printAssignment = [&]()
	{
	std::cout << "supertile " << supertileId << ": ";
	for(const auto e : (*dbg_supertileAssignments)[supertileId])
	  std::cout << "<" << e << ">";
	std::cout <<  std::endl;
	};
	std::cout << "BEFORE ";
	printAssignment();
	
	if(add)
	  {
	    const bool success = (*dbg_supertileAssignments)[supertileId].insert(tileId).second;
	    assert(success);
	  }
	else
	  {
	    const auto cnt = (*dbg_supertileAssignments)[supertileId].erase(tileId);
	    hassertm2( cnt==1, supertileId, tileId);
	  }
#endif
	
	const auto oldCnt = supertilesCnts[supertileId];
	hassertm(add || oldCnt>0, supertileId);
	const int mod=(add ? 1 : -1);
	const auto newCnt = oldCnt+mod;	  

	auto sv = &supertiles[supertileId*nElems_tile];
	
	for(uint32_t i=0; i<nElems_tile; i++)
	  {
	    auto & s = sv[i];

	    s = getUpdatedValue(s, oldCnt, mod, tiles[tileId*nElems_tile+i]);
	    hassertm(sanityCheckValue(s), s);
	    // if(newCnt==0)
	    //   s*=0.;
	    // else
	    //   s=(oldCnt*s+mod*tiles[tileId*nElems_tile+i])/newCnt;
	  }	

	supertilesCnts[supertileId]=newCnt;

	#ifdef DBG_SUPERTILES_DF
	std::cout << "supertile cnt of " << supertileId << " changed from " << oldCnt << " to " << newCnt << std::endl;

	std::cout << "AFTER ";
	printAssignment();
	#endif
      }
      
    public:

      template<typename E>
#ifdef __NVCC__
      __host__ __device__
#endif
      auto target2supertileId(E targetId) const
      {
	hassertm3(targetId >= nElems_source
		  && targetId-nElems_source<nElems_target,
		  targetId, nElems_source, nElems_target);
	return targetId-nElems_source;
      }
	
      
      template<typename E>
#ifdef __NVCC__
      __host__ __device__
#endif
      void swapSequence(E* swapTiles, unsigned int n)
      {
	
	assert(tiles != 0);
	assert(supertiles != 0);

	auto sanityCheck = [this](E e)
	{
	  hassertm2(e.x < nElems_source, e.x, nElems_source);
	  hassertm3(e.y >= nElems_source
		    && e.y-nElems_source<nElems_target,
		    e.y, nElems_source, nElems_target);
	};

	//std::cout << "|"<< nElems_source << " " << nElems_target << std::endl;
	  
	//std::cout << "---------------------- " << n << "\n";

	for(unsigned int i=0; i<n; i++)
	  {
	    const E target=swapTiles[i];
	    
	    const E source=swapTiles[(i+1)==n ? 0 : (i+1)];
	    
	    sanityCheck(source);
	    sanityCheck(target);

#ifdef DBG_SUPERTILES_DF
	    std::cout << "assignment of tile <" << source.x << "> is changed from supertile " << source.y << "(" << target2supertileId(source.y) << ") to supertile " << target.y << "(" << target2supertileId(target.y) << ")" << std::endl;
#endif	    

	    auto it_from = qt(target2supertileId(source.y));
	    auto it_to = qt(target2supertileId(target.y));

	    // only do it on the lowest level for debugging purposes
	    updateQTLevel(it_from.idx, source.x, false);
	    updateQTLevel(it_to.idx, source.x, true);
	    
	    while(true)
	      {
		//advance first, we're not interested in leaves		
		auto success0 = it_from();
		auto success1 = it_to();
		assert(success0 == success1);

		// no update if we reached root
		// OR we reached the point where the
		// swap does not change anything on a certain level
		if(!success0 || it_from.idx == it_to.idx)
		  break;
		 
		// source loses the tile, target gets new tile
		updateQTLevel(it_from.idx, source.x, false);
		updateQTLevel(it_to.idx, source.x, true);
	      }
	  }
	
      }


      template<typename A>
      auto createMap_assignments_host(A* assignments)
      {
	std::cout << "begin " << __func__ << std::endl;
	std::vector<unsigned int> aCnts(nElems_target, 0);

	std::vector<uint32_t> tileAssignment(assignments->size());
	  
	for(size_t i=0; i<assignments->size(); i++)
	  {
	    assert((assignments->cend(i)-assignments->cbegin(i))==1);
	    std::cout << i;
	    for(auto it = assignments->cbegin(i);
		it != assignments->cend(i); it++)
	      {		
		auto supertileId = assignments->getIdx(*it);
		assert(supertileId >= nElems_source);
		supertileId -= nElems_source;
		aCnts[supertileId]++;
		tileAssignment[i] = supertileId;
		std::cout << " [" << supertileId << ":" << assignments->getMass(*it) << "]";
	      }
	    std::cout << std::endl;	    
	  }

	for(size_t i=0; i<nElems_target;i++)
	  hassertm3(supertilesCnts[i]==aCnts[i],i,
		    supertilesCnts[i],aCnts[i]);
	std::cout << "end " << __func__ << std::endl;
	return tileAssignment;
      }

      V2<uint32_t> tileGridDim;
      const TILE* tiles;
      TILE* supertiles=0;
      unsigned int* supertilesCnts=0;
      uint32_t nElems_tile=0;

      const unsigned int* qtNeighbors;
#ifdef DBG_SUPERTILES_DF
      std::vector<std::set<uint32_t>>* dbg_supertileAssignments=0;
#endif
    };

#endif // __INTERVOL_DIST_FUNC__
