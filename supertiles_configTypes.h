#ifndef __SUPERTILES_CONFIG_TYPES__
#define __SUPERTILES_CONFIG_TYPES__

enum distFuncType_t
  {
    distFuncType_norm2=1,
    distFuncType_cosine_normalized=2,
    distFuncType_none=0
  };

// enum filterType_t
//   {
//     filterType_box=1,
//     filterType_gaussian_elem=2,
//     filterType_none=0
//   };

enum distVisType_t
  {
    distVisType_separator=1,
    distVisType_none=0,
    distVisType_posScale=2
  };

enum repAggregationType_t
  {
    repAggregationType_front=0,
    repAggregationType_average=1,
    repAggregationType_central=2,
    repAggregationType_mcmc=3,
    repAggregationType_chart=4,
    repAggregationType_average_colRGB=5,
  };


#endif //__SUPERTILES_CONFIG_TYPES__
