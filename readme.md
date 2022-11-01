<div id="top"></div>


<!-- ABOUT THE PROJECT -->
## About The Project

This is the implementation of grid layout method described in this EuroVis 2022 paper:

S. Frey, “Optimizing Grid Layouts for Level-of-Detail Exploration of Large Data Collections,” 2022, <https://doi.org/10.1111/cgf.14537>.

It has been tested on macOS BigSur and Ubuntu 22.04. 


## Prerequisites

These libraries are prerequisites for building

* Cairo graphics library

* bzip2

* cmake


<!-- GETTING STARTED -->
## Getting Started

*ldg_core* can simply be built via cmake.
These are exemplary steps to get you started on Linux/macOS.

* Obtain the code from the repository https://github.com/freysn/ldg_core , e.g., via <code>git clone https://github.com/freysn/ldg_core.git</code>

* Enter the directory: <code>cd ldg_core</code>

* Create a build directory: <code>mkdir build</code>

* Now create a Makefile using cmake <code>cmake ..</code>. Some variants of this might be useful:

- cmake -DCMAKE_BUILD_TYPE=Release ..
- cmake -DCMAKE_BUILD_TYPE=Debug ..
- cmake .. -DCMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=clang-mp-14 -D CMAKE_EXE_LINKER_FLAGS=-lc++ (Note that on macOS, the standard compiler can be used but offers no support for OpenMP. For OpenMP, a different compiler needs to be used. For example, on my mac, I use clang-14 (obtained via macports).).

<!-- USAGE EXAMPLES -->
## Usage

*ldg_core* can both be used to optimize the placement of member of data collections as well as rendering the resulting grid.

- run_examples_1k.sh exemplifies the usage for datasets with 1024 members
- run_examples_full.sh demonstrates the usage for the full datasets (note that this requires more compute time and memory to run)

The required data package can be downloaded from here: <https://www.dropbox.com/s/cp4rsq8cyjfl7qs/ldg_data.zip?dl=0>.

Unless any termination criteria are specified, the optimizer runs until it is interrupted by CTRL-C.  
By default, it prints a short update with the current evaluation score and for how long the assignment has been unchanged (among others). 
It then completes the ongoing pass phase and then writes the assignment file "qtLeafAssigning.raw.bz" in the specified output directory. 

<!-- Please see the code snippet for an example how an assignment file asFile can be read and a resulting image can be generated (by means of rearranged RGB colors in an image img with image dimensions imgDim): -->

<!-- ``` -->
<!-- const auto qtLeafAssignment= -->
<!-- 	supertiles::place::read_qtLeafAssignment(asFile, -->
<!-- 						 imgDim); -->

<!-- const auto leaf2gridPos = -->
<!-- 	helper::invertMap(gridPos2QTLeaves(imgDim)); -->


<!-- std::vector<unsigned char> o(helper::ii2n(imgDim)*nChannels, 0); -->
<!-- for(const auto leafPos : helper::range_n(helper::ii2n(imgDim))) -->
<!-- 	{ -->
<!-- 	  const auto imgIdx = qtLeafAssignment[leafPos]; -->
<!-- 	  const auto gridPos = leaf2gridPos[leafPos]; -->
<!-- 	  for(const auto c : helper::range_n(nChannels)) -->
<!-- 	    o[nChannels*gridPos+c] = img[nChannels*imgIdx+c]; -->
<!-- 	} -->
      
<!-- helper::cimgWrite("colRGB_fromPNG.png", &o[0], imgDim, nChannels); -->
<!-- ``` -->


### command line switches
Useful command line switches for ldg_core (incomplete list, but should contain most important ones to get started); they all expect one argument:

* --dataconf, -d | custom config format describing the input dataset (see description below).
* --outDir | directory where to store generated assignment files
* --distFuncType | 1: Euclidean distance or 2: Cosine distance
* --seed | seed for the pseudo RNG
* --nNeighbors | 4 (for 4-neighborhood, default) or 8 (for 8-neighborhood)
* --termIterations | specify how many optimization iterations the code should run (default: infinity)
* —termSame” | interrupt after the assignment hasn't been changed for said number of iterations.
* --termTime | interrupt after running for a certain number of seconds


### data config files

_config_ files contain textual descriptions of the files to load. It assumes that the data is stored in a raw file.

* VOLUME\_FILE _arg_ | name of the raw file to load
* VOLUME\_DIM _arg0_ _arg1_ _arg2_ | volume dimension in x, y, and z direction (integers)
* VOLUME\_DATA_TYPE _arg_ | element type, possible values: _UCHAR_, _USHORT_, _FLOAT_, _DOUBLE_, _UCHAR4_
* VOLUME\_FORMAT | data format, possible values: _RAW_, _RAW_BZ2_ (raw data file compressed with bzip2)

<!-- _For more examples, please refer to the [Documentation](https://example.com)_ -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- ROADMAP -->
<!-- ## Roadmap -->

<!-- - [] Feature 1 -->
<!-- - [] Feature 2 -->
<!-- - [] Feature 3 -->
<!--     - [] Nested Feature -->

<!-- See the [open issues](https://github.com/github_username/repo_name/issues) for a full list of proposed features (and known issues). -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- CONTRIBUTING -->
<!-- ## Contributing -->

<!-- Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**. -->

<!-- If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement". -->
<!-- Don't forget to give the project a star! Thanks again! -->

<!-- 1. Fork the Project -->
<!-- 2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`) -->
<!-- 3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`) -->
<!-- 4. Push to the Branch (`git push origin feature/AmazingFeature`) -->
<!-- 5. Open a Pull Request -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



## License

Distributed under the MIT License. See `license.md` for more information.

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- CONTACT -->
## Contact

Steffen Frey - s.d.frey@rug.nl

Project Link: [https://github.com/freysn/ldg_core](https://github.com/freysn/ldg_core)

An interactive demo for exploring LDGs is provided here: <https://ldg-demo.github.io>.

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- ACKNOWLEDGMENTS -->
<!-- ## Acknowledgments -->

<!-- * []() -->
<!-- * []() -->
<!-- * []() -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
<!-- [contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge -->
<!-- [contributors-url]: https://github.com/github_username/repo_name/graphs/contributors -->
<!-- [forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge -->
<!-- [forks-url]: https://github.com/github_username/repo_name/network/members -->
<!-- [stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge -->
<!-- [stars-url]: https://github.com/github_username/repo_name/stargazers -->
<!-- [issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge -->
<!-- [issues-url]: https://github.com/github_username/repo_name/issues -->
<!-- [license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge -->
<!-- [license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt -->
<!-- [linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555 -->
<!-- [linkedin-url]: https://linkedin.com/in/linkedin_username -->
<!-- [product-screenshot]: images/screenshot.png -->
