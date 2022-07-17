#include "supertiles_place.h"

// bool isNumber(const std::string& s)
// {
//     return !s.empty() && std::find_if(s.begin(), 
//         s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
// }

int main(int argc, const char** argv)
{
  // std::cout << helper::apprEq(static_cast<double>(0.104782), static_cast<double>(0.1047821)) << std::endl;;
  
  // return 0;
#ifndef NO_OMP
  if(false)
    omp_set_num_threads(1);
#endif
  
  // const auto po=initPO(argc, argv);

  // std::string cfg;

  // if(!po.configFNames.empty())
  //   cfg=po.configFNames.front();
  // else
  //   cfg="65536";

  // std::cout << "CFG: " << cfg << std::endl;
  // if(isNumber(cfg))
  //   {
  //     std::cout << "start benchmark mode\n";
  //     supertiles::place::benchmark_rgb(po, std::stoi(cfg));
  //     //supertiles::place::benchmark_scalar();
  //     //supertiles::place::testSupertilesIO_rgb();
  //   }
  // else
  const std::string e(argv[0]);
  if(e.size()>=3 && std::string("f32")==e.substr(e.size()-3, 3))
    {
      std::cout << "f32" << std::endl;
      supertiles::place::run_f32(argc, argv);
    }
  else
    {
      std::cout << "f64" << std::endl;
      supertiles::place::run(argc, argv);
    }
  return 0;
}

