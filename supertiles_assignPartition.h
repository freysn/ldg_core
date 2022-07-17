// resume process sequence until sequence is processed to the end

// const unsigned int mode_doExchange = 1;
// const unsigned int mode_doExchangeLog = 2;


// TODO: new variant that solves full assignment problem (e.g, with hungarian algorithm) for provided subset; for this, always subtract all tiles in plan from supertiles; IDEA: COULD ACTUALLY DO THE SUBTRACTION BEFOREHAND, THEN JUST ADD NEW STUFF AFTER CHECKS ARE COMPLETE (not very efficient, but simple)

template<unsigned int mode, typename F, typename A, typename PLAN, typename DIST_FUNC,/* typename P,*/ typename IDX, typename CACHE>
#ifdef USE_NVCC
    __device__
#endif
  thrust::tuple<F,uint32_t, uint32_t> supertiles_assignPartition(A& assignments, PLAN& plan, DIST_FUNC distFunc, /*const P& pos,*/ const IDX& lidx, CACHE v, uint32_t planElemOff=0,
						      thrust::tuple
						      <typename A::e_t*,unsigned int*, unsigned int>
						      exchangeLog
						      =thrust::make_tuple((typename A::e_t*)0,
									  (unsigned int*)0, 0))
  {
    //std::cout << "-----\n";
#ifndef USE_NO_BVOL
    typedef typename std::remove_reference<decltype(plan._bvols[0])>::type BVOL;
#endif
    typedef typename std::remove_reference<decltype(plan.elems[0])>::type E;
    typedef typename std::remove_reference<decltype(assignments.begin(0))>::type A_it;
    typedef typename std::remove_const
      //<typename std::remove_reference<decltype(distFunc.pos_source[0])>::type>::type F3;
      <typename std::remove_reference<decltype(distFunc.getPosSource(0))>::type>::type F3;

    F d_out = 0.;

    uint32_t mass =
      //std::numeric_limits<uint32_t>::max()
      UINT32_MAX
      ;
    
    IDX nElemsListMax;
    if(!plan.getNElemsList(nElemsListMax, lidx))
      return thrust::make_tuple(0., 0, 0);

    E e_back;

    if(!plan.get(e_back, lidx, planElemOff))
      {
        assert(false);
        return thrust::make_tuple(0., 0, 0);
      }    

    //F3 p_back = distFunc.pos_source[e_back];
    F3 p_back = distFunc.getPosSource(e_back);

    E e_front = e_back;
    F3 p_front = p_back;
    
    //
    // prev best is what is taken from the predecessor source element
    //
    A_it prevBest = static_cast<A_it>(0);
    
    IDX nElemsList=0;

#ifndef USE_NO_BVOL
    BVOL bvol_back = plan._bvols[e_back];
#endif
    
    {
      //when re-entering, we do not resume with old i (nElemsOff-1) because exchange information needs updating
    uint32_t i=planElemOff+1;
    while(true)
      {
	if(i==nElemsListMax)
	  return thrust::make_tuple(0., 0, 0);
	assert(i<nElemsListMax);
        E e_backm1 = e_back;
        
        if(!plan.get(e_back, lidx, i))
	  {
	    // this should not happen I think
	    assert(false);
	    return thrust::make_tuple(0., 0, 0);
	  }

        //std::cout << "exchange " << i << ": " << e_back << std::endl;

        const F3 p_backm1 = p_back;
        //p_back = distFunc.pos_source[e_back];
	p_back = distFunc.getPosSource(e_back);

	//
	// BVOL TEST END
	//
	bool exchangePossible = true;
	bool continueSequencePossible = true;

#ifndef USE_NO_BVOL
	BVOL bvol_backm1 = bvol_back;
	bvol_back = plan._bvols[e_back];
	
	const F d_out_bestCase_one_delta
	  = potentialImprovement<F>(p_backm1, p_back, bvol_backm1);
	const F d_out_bestCase_two_delta
	  = potentialImprovement<F>(p_back, p_front, bvol_back);

	continueSequencePossible = (d_out + d_out_bestCase_one_delta < 0.);
	exchangePossible = (d_out + d_out_bestCase_one_delta + d_out_bestCase_two_delta < 0.);
#endif

	//exchangePossible = true;
	
	//
	// BVOL TEST BEGIN
	//

	if(continueSequencePossible || exchangePossible)
	  {
#ifndef USE_NO_BVOL
	    BVOL bvol;
	    bvol.prepareUpdate();
#endif
	    //std::cout << "current assignments " << e_backm1 << ": ";

	    //
	    // identify best target element to transfer from p_backm1 to p_back
	    //
	    {
	      FindBest_empty_t<A_it, F> rslt;
	      {	  
		A_it abegin = assignments.begin(e_backm1);
		//continue;
		//A_it aend = assignments.end(e_backm1);
		A_it aend = abegin + assignments.size(e_backm1);
	  
		rslt =
		  INTERVOL_NS::findBest_empty<A,F>(abegin,
						 aend,
#ifdef SUPERTILES
						   assignments.begin(e_back),				
#endif
						 p_backm1,
						 // give to
						 p_back,
						 distFunc,
						 //pos,						 
#ifndef USE_NO_BVOL
						 bvol,
#endif
						 prevBest);
	      }
	
	      //std::cout << std::endl;
#ifndef USE_NO_BVOL
	      if(bvol_backm1 != bvol)
		{
		  plan._bvols[e_backm1] = bvol;
		  //bvol_backm1 = bvol;
		}
#endif
	      //
	      // could not find a target element in p_backm1 for p_back
	      // I think this happens when it just has one type of target element
	      // and it gets exactly this target as input (?)
	      //
	      if(rslt.best == static_cast<A_it>(0))
		{
		  return thrust::make_tuple(0., 0, 0);
		}
        
	      d_out += rslt.best_diff;

	      //std::cout << "BVOLTEST " << d_out << " vs " << d_out_bestCase_one << std::endl;
#ifndef USE_NO_BVOL
	      /*
#ifndef NDEBUG
	      if(rslt.best_diff+0.00001 <= d_out_bestCase_one_delta)
		printf("error %f %f\n", rslt.best_diff, d_out_bestCase_one_delta);
#endif
	      */
	      assert(rslt.best_diff+0.00001 >= d_out_bestCase_one_delta);
#endif

	      prevBest = rslt.best;

	      /*
		if(planElemOff > 0 && i>planElemOff)
		assert(v[i].best_diff == rslt.best_diff);
	      */

	      assert(i< plan.getMaxNElemsList());
	      v[i] = rslt;

	
	
	      mass =
#ifndef USE_NVCC
		std::
#endif
		min(mass, A::getMass(*rslt.best));
	    }
	  }
	else
	  d_out = 1.;
	
	assert(mass > 0);
	{

#ifndef USE_NO_BVOL
	  if(exchangePossible && (d_out + d_out_bestCase_two_delta < 0.))
#endif
	    {
	      //
	      // identify best target element to transfer from p_back to p_front
	      //
	      A_it abegin = assignments.begin(e_back);
	      //continue;
	      //A_it aend = assignments.end(e_backm1);
	      A_it aend = abegin + assignments.size(e_back);

#ifdef SUPERTILES
	      assert(assignments.size(e_back)==1);
#endif

#ifndef USE_NO_BVOL
	      BVOL bvol;
	      bvol.prepareUpdate();
#endif
	  
	      FindBest_empty_t<A_it, F> rslt =
		INTERVOL_NS::findBest_empty<A,F>(abegin,
					       aend,
#ifdef SUPERTILES
						 assignments.begin(e_front),				
#endif
					       p_back,
					       // give to
					       p_front,
					       distFunc,
					       //pos,
#ifndef USE_NO_BVOL
					       bvol,
#endif
					       prevBest);

#ifndef USE_NO_BVOL
	      if(bvol != bvol_back)
		{
		  plan._bvols[e_back] = bvol;
		  bvol_back = bvol;
		}
#endif

	      //
	      // check if valid sequence has been found
	      //
	      if(rslt.best != static_cast<A_it>(0))
		{
		  const auto d_final = d_out + rslt.best_diff;

#ifndef USE_NO_BVOL		  
		  assert(rslt.best_diff+0.00001 >= d_out_bestCase_two_delta);
#endif
		  if(d_final < 0.)
		    {
		      assert(planElemOff< plan.getMaxNElemsList());
		      v[planElemOff] = rslt;
		      nElemsList = i+1;
		      d_out = d_final;
		      assert(mass > 0);

		      const auto curMass = A::getMass((*rslt.best));
		      assert(curMass > 0);
		      mass =
#ifndef USE_NVCC
			std::
#endif
			min(mass, curMass);

		      assert(mass > 0);
		      prevBest = rslt.best;

		      //printf("%d %d: %d %d\n", threadIdx.x, lidx, e_front, e_back);
		      //std::cout << __PRETTY_FUNCTION__ << " found new sequence, l=" << nElemsList-planElemOff << " " << d_out << " " << mass << " " << planElemOff << " " << nElemsListMax << std::endl;
		      break;
		    }	      
		}
	    }
	  //
	  // check if sequence is not promising (i.e., I do not gain from exchanges so far)
	  //
	  if(d_out > 0.)
	    {
	      //
	      // abort sequence
	      //
	      prevBest = static_cast<A_it>(0);
	      d_out = 0.;
	      planElemOff = i;
	      mass = UINT32_MAX;
	      p_front = p_back;
	      e_front = e_back;
	    }
	}
	i++;

	
        //std::cout << "best candidate for " << e_back << " from " << e_backm1 << ": "<< getIdx(*rslt.best) << " " << getMass(*rslt.best) << std::endl;
      }
    }
    assert(mass>0);    
    
    //
    // swap does not improve assignment quality
    //
    //if(d_out >= 0.)
    //return 0.;
    assert(d_out < 0.);

    //return 666.;

    //std::cout << "swap with " << nElemsList << " elements " << d_out << std::endl;;
    /*
    if(nElemsList > 2)
      {
        std::cout << "swap with " << nElemsList << " elements\n";
        //assert(false);
      }
    */
    //
    // CHECK IF ZEROTH ELEMENT ALREADY HAS
    // MASS WITH ID THAT IS ASSIGNED TO IT
    //    
    {
      //E e;
      //plan.get(e, lidx, nElemsList-1);
      assert(v[planElemOff+1].other_match_it == 0);
      v[planElemOff+1].other_match_it =
        INTERVOL_NS::findFirstId<A>(assignments.begin(e_front),
                              assignments.end(e_front),
				 A::getIdx(*(v[planElemOff].best)));
    }    

    
    unsigned int exchangeLogId = 0;

    auto exchangeLogAdd = [&exchangeLog, &exchangeLogId](unsigned int alpha, unsigned int omega_received)
      {	
	//uint2 info;
	typename A::e_t info;
	A::setIdx(info, alpha);
	A::setMass(info, omega_received);
	//std::cout << "exchangeLog_add " << exchangeLogId << " " << info.x << " " << info.y << std::endl;
	assert(exchangeLogId < *thrust::get<1>(exchangeLog));
	thrust::get<0>(exchangeLog)[exchangeLogId] = info;
	exchangeLogId++;
      };
      
    if(mode & mode_doExchangeLog)
      {
	const auto nElemsSequence = nElemsList-planElemOff;
	assert(nElemsSequence > 1);
	assert(thrust::get<1>(exchangeLog) != 0);
#if defined(USE_CUDA_RUN_DEVICE)
	exchangeLogId = atomicAdd(thrust::get<1>(exchangeLog), nElemsSequence+1);
#else
#ifdef USE_OMP_INTERVOL
	#pragma omp atomic capture
#endif
	{
	  exchangeLogId = *thrust::get<1>(exchangeLog);
	  *thrust::get<1>(exchangeLog) += (nElemsSequence+1);
	}
#endif
	//std::cout << "-----------------\n";
	exchangeLogAdd(nElemsSequence, mass);
	
	assert(*thrust::get<1>(exchangeLog) <= thrust::get<2>(exchangeLog));
      }
    
    //
    // APPLY SWAPS
    //
    auto prev_best = A::getIdx(*v[nElemsList-1].best);

#ifdef SUPERTILES
    std::vector<uint3> swapElems;
    //std::cout << "ASSIGNMENTS BEFORE SWAP " << assignments << std::endl;
#endif
    {
      auto k_prev = nElemsList-1;
      for(IDX k=planElemOff; k<nElemsList; k++)
	{
	  assert(v[k].best != static_cast<A_it>(0));
	  typename A::e_t n;

	  //const auto k_prev = (k+nElemsList-1)%nElemsList;
	  A::setIdx(n, prev_best);

	  prev_best = A::getIdx(*v[k].best);
          
	  A::setMass(n, mass);
        
	  E e;
	  plan.get(e, lidx, k_prev);
#ifdef SUPERTILES
	  swapElems.push_back(make_uint3(e, v[k].best->x, v[k].best->y));
#endif
	  
	  
	  //std::cout << "e: " << e << " " << v[k].best->x << " " <<  v[k].best->y << std::endl;
	    
        
	  assert(v[k].other_match_it == static_cast<A_it>(0)
		 || A::getIdx(*v[k].other_match_it) == A::getIdx(n));

	  if(mode & mode_doExchange)
	    {
	      INTERVOL_NS::applySwap(assignments,
				     e, v[k],
				     v[k].other_match_it, n);


#ifndef USE_NO_BVOL
	      const auto targetIdx = A::getIdx(n);

	      auto rslt = distFunc.getPosTarget_inRange(targetIdx);
	      assert(thrust::get<1>(rslt));	    
	      plan._bvols[e].update(thrust::get<0>(rslt));
	      
	      //const auto localTargetIdx = distFunc.localTargetIdx(targetIdx);
	      //plan._bvols[e].update(distFunc.pos_target[localTargetIdx]);
#endif
	    }

	  if(mode & mode_doExchangeLog)
	    exchangeLogAdd(e, A::getIdx(n));
	  
	  k_prev = k;
	}
    }

#ifdef SUPERTILES
    distFunc.swapSequence(&swapElems[0], swapElems.size());
#endif
    //std::cout << /*__PRETTY_FUNCTION__ <<*/ " d_out "  << d_out << " mass " << mass << " " << d_out*mass << std::endl;
    //assert(mass==1);
    
    return thrust::make_tuple(d_out*mass, planElemOff, nElemsList);
  }
