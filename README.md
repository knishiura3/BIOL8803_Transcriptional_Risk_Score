<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

![Anaconda](https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white)
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">

<h3 align="center">GWAS & eQTL Colocalization</h3>

  <p align="center">
    project_description
    <br />
    <a href="https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score">View Demo</a>
    ·
    <a href="https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/issues">Report Bug</a>
    ·
    <a href="https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project



<!-- Here's a blank template to get started: To avoid retyping too much info. Do a search and replace with your text editor for the following: `knishiura3`, `BIOL8803_Transcriptional_Risk_Score`, `twitter_handle`, `linkedin_username`, `email_client`, `email`, `project_title`, `project_description` -->

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ### Built With

* [![Next][Next.js]][Next-url]
* [![React][React.js]][React-url]
* [![Vue][Vue.js]][Vue-url]
* [![Angular][Angular.io]][Angular-url]
* [![Svelte][Svelte.dev]][Svelte-url]
* [![Laravel][Laravel.com]][Laravel-url]
* [![Bootstrap][Bootstrap.com]][Bootstrap-url]
* [![JQuery][JQuery.com]][JQuery-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p> -->



<!-- GETTING STARTED -->
## Getting Started
<!-- 
This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps. -->

### Prerequisites

<!-- This is an example of how to list things you need to use the software and how to install them. -->
* anaconda/mamba (only tested on Ubuntu 22.04.1 LTS in WSL2. Things might be different for Windows/Mac/non-virtual Linux)

### Installation


1. Clone the repo
   ```
   git clone https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score
   ```
2. Create a conda environment using the yaml file provided
    ```
    conda env create -f=/path/to/<env yml file>
    ```
2. Open R within the conda environment & install R packages (includes packages needed for R kernel usage in jupyter nb)

* R 
  ```
  install.packages(c("devtools", "shiny", "shinythemes", "shinycssloaders", "DT", "slickR", 
                 "duckdb", "fs", "tidyverse", "DBI", "glue", "dplyr", "coloc", "ggplot2", "httpgd"))
  devtools::install_github("jrs95/gassocplot")
  devtools::install_github("mrcieu/gwasglue")
  
  install.packages(c("devtools",
    "languageserver",
    "BiocManager",
    "tidyverse",
    "tictoc",
    "httpgd"), dependencies = TRUE)
  BiocManager::install(version = "3.14") 
  devtools::install_github("IRkernel/IRkernel", dependencies = TRUE)
  IRkernel::installspec()

  ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage
1) Build the parquet from eQTLGen Summary statistics.
2) Query the parquet using the desired GWAS ID from openGWAS.
Refer to https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/blob/parquet_duckDB/Coloc_tutorial.ipynb for example usage

A shiny app implementation is also available at https://genapp2022.biosci.gatech.edu/team1/

<!-- Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources. -->

<!-- _For more examples, please refer to the [Documentation](https://example.com)_ -->

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Kenji Nishiura - kenji@gatech.edu

Colin Naughton - Naughtoncolin@gmail.com

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* Many thanks to [Kenji Gerhardt](https://github.com/KGerhardt) for his assistance with coding best practices and optimization

## Team Members

* Andy Chea
* Colin Naughton
* Kenji Nishiura
* Jasmyn Pellebon

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/knishiura3/BIOL8803_Transcriptional_Risk_Score.svg?style=for-the-badge
[contributors-url]: https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/knishiura3/BIOL8803_Transcriptional_Risk_Score.svg?style=for-the-badge
[forks-url]: https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/network/members
[stars-shield]: https://img.shields.io/github/stars/knishiura3/BIOL8803_Transcriptional_Risk_Score.svg?style=for-the-badge
[stars-url]: https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/stargazers
[issues-shield]: https://img.shields.io/github/issues/knishiura3/BIOL8803_Transcriptional_Risk_Score.svg?style=for-the-badge
[issues-url]: https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/issues
[license-shield]: https://img.shields.io/github/license/knishiura3/BIOL8803_Transcriptional_Risk_Score.svg?style=for-the-badge
[license-url]: https://github.com/knishiura3/BIOL8803_Transcriptional_Risk_Score/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/kenji-nishiura/
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
