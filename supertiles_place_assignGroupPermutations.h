#ifndef __SUPERTILES_PLACE_ASSIGN_GROUP_PERMUTATIONS__
#define __SUPERTILES_PLACE_ASSIGN_GROUP_PERMUTATIONS__

#include "supertiles_place_cost.h"
#include "supertiles_place_init.h"
#include <unordered_map>

#include "helper_hungarian.h"

#define ASSIGN_GROUP_HUNGARIAN
//#define ASSIGN_GROUP_PERMUTATIONS

namespace supertiles
{
  namespace place
  {

    template<typename T>
    T factorial(T n)
    { 
      // single line to find factorial 
      return (n==1 || n==0) ? 1: n * factorial(n - 1);  
    }
    
    auto createAssignPermutations(size_t groupSizeMax)
    {
      using idx_t = uint32_t;
      using idcs_t = std::basic_string<idx_t>;

      std::vector<std::vector<idcs_t>> assignmentCandidatesv(groupSizeMax+1);


      for(const auto & N : helper::range_be(2, groupSizeMax+1))
	{

	  std::vector<idcs_t> assignmentCandidates;
   
   
	  {
	    idcs_t assignment;
   
	    for(uint32_t i=0; i<N; i++)
	      {
		assignment.push_back(i);
	      }
	    //std::cout << std::endl;

     

	    assignmentCandidates.reserve(factorial(N));

	    do {
	      assignmentCandidates.push_back(assignment);
	    }
	    while(std::next_permutation(assignment.begin(), assignment.end()));

	    assert(assignmentCandidates.size()==factorial(N));
	  }

	  assignmentCandidatesv[N]=assignmentCandidates;
	}
      return assignmentCandidatesv;
    }

#ifdef __NO_OMP
#error "SHOULD NOT BE DEFINED"
#endif     

#ifdef ASSIGN_GROUP_HUNGARIAN
#define __NO_OMP
    //#include "supertiles_place_assignGroupPermutations_.h"
#include "supertiles_place_assignGroupHungarian_.h"
#undef __NO_OMP
#endif

#ifdef ASSIGN_GROUP_PERMUTATIONS
namespace omp
{
#include "supertiles_place_assignGroupPermutations_.h"
}
#endif

  }
}

#endif //__SUPERTILES_PLACE_ASSIGN_GROUP_PERMUTATIONS__
