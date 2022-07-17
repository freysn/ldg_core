<div id="top"></div>


<!-- ABOUT THE PROJECT -->
## About The Project

This is the implementation of grid layout method described in this paper EuroVis 2022 paper:

S. Frey, “Optimizing Grid Layouts for Level-of-Detail Exploration of Large Data Collections,” 2022, doi: 10.1111/cgf.14537.

It has been tested on macOS BigSur and Ubuntu 20.04. 
<!-- Below, a quick sketch of the main concept is provided. Please refer to the paper for details. -->

<!-- * TODO -->



<!-- ### Built With -->

<!-- * [Next.js](https://nextjs.org/) -->
<!-- * [React.js](https://reactjs.org/) -->
<!-- * [Vue.js](https://vuejs.org/) -->
<!-- * [Angular](https://angular.io/) -->
<!-- * [Svelte](https://svelte.dev/) -->
<!-- * [Laravel](https://laravel.com) -->
<!-- * [Bootstrap](https://getbootstrap.com) -->
<!-- * [JQuery](https://jquery.com) -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- GETTING STARTED -->
<!-- ## Getting Started -->

<!-- <\!-- This is an example of how you may give instructions on setting up your project locally. -\-> -->
<!-- To get a local copy up and running follow these simple example steps. -->

## Prerequisites

These libraries are prerequisites for building

* Steffen's helper library
```git clone git@github.com:freysn/helper.git
```

* Cairo graphics library (could also be removed if needed)

* bzip2

<!-- ### Installation -->

<!-- 1. Get a free API Key at [https://example.com](https://example.com) -->
<!-- 2. Clone the repo -->
<!--    ```sh -->
<!--    git clone https://github.com/github_username/repo_name.git -->
<!--    ``` -->
<!-- 3. Install NPM packages -->
<!--    ```sh -->
<!--    npm install -->
<!--    ``` -->
<!-- 4. Enter your API in `config.js` -->
<!--    ```js -->
<!--    const API_KEY = 'ENTER YOUR API'; -->
<!--    ``` -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- USAGE EXAMPLES -->
## Usage

### command line switches
Useful command line switches to _timeStepSelector_ : 

* --dataconf, -d _arg_ (default: testData/bottle.config) | custom config format describing the input dataset (see description below).

* --maxSelect, -n _arg_ (default: 16) | maximum number of time steps to be selected (the creates selections for [1,2,3, .., maxSelect] time steps).

* --batch,-b _arg_ (default: "selections") | output file in which the selections are stored (one line per number of time steps, indices separated by white space).

* --distCache _args_ (default: "") | cache file for computed distances. If it exists, distances are loaded from this file, otherwise they are newly computed and stored in the specified file.

* --selectionFile _arg_ (default: none) | expects the file name of a one-line text file with whitespace-separated indices of time steps. Outputs a file of the same name with postifx \_eval (e.g., "selection" results in "selection\_eval"). _Note: if provided, no selections are computed, but only the given selection is evaluated._

* --fixFirstLastTimeStep (default: false) | Forces selection of the first and last time step.

### data config files

_config_ files contain textual descriptions of the files to load. It assumes that each time step is stored in an individual file (see _testData/bottle.conf_ for an example).

* VOLUME\_FILE _arg_ | names of the volume files, wildcard with * supported
* VOLUME\_DIM _arg0_ _arg1_ _arg2_ | volume dimension in x, y, and z direction (integers)
* VOLUME\_DATA_TYPE _arg_ | element type, possible values: _UCHAR_, _USHORT_, _FLOAT_, _DOUBLE_
* VOLUME\_FORMAT | data format, possible values: _PNG_, _RAW_, _RAW_BZ2_ (raw data file compressed with bzip2)
* VOLUME\_NUM_TIMESTEPS _arg_ | number of time steps to load

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



<!-- LICENSE -->
<!-- ## License -->

<!-- Distributed under the MIT License. See `LICENSE.txt` for more information. -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- CONTACT -->
<!-- ## Contact -->

<!-- Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email@email_client.com -->

<!-- Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name) -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- ACKNOWLEDGMENTS -->
<!-- ## Acknowledgments -->

<!-- * []() -->
<!-- * []() -->
<!-- * []() -->

<!-- <p align="right">(<a href="#top">back to top</a>)</p> -->



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
