# 1<sup>ST</sup> ORDER MRI SIMS
## First-order Spatial Encoding Simulations for Improved Accuracy in the Presence of Strong $\mathrm{B}_0$ and Gradient Field Variations

[[`Paper`]()] [[`Github`](XXLINK)] [[`BibTeX`](#cite)]


## Overview

This repository contains code to perform 1<sup>st</sup>-order spatial encoding MRI simulations, developed especially for highly variable and inhomogeneous main magnetic (B<sub>0</sub>) and gradient fields (BG<sub>xyz</sub>). This simulation code is a starting point for exploring how to capture the encoding effects of strong field variations accurately and efficiently to enable the assessment of distortions, signal dropout, and foldover artifacts.

The following figure determines an example scenario in which these simulations can be used for improved accuracy:

![why1st][https://github.com/radhikatibrewala/1stOrderMRISims/blob/master/image.png]

For more details, please refer to the paper.

## Installation

The code requires `julia >= v"1.9.0"`

Install 1<sup>ST</sup> ORDER MRI SIMS: clone the repository locally and install with

```
git clone https://github.com/radhikatibrewala/1stOrderMRISims.git
```

## Usage
The repository is centered around the ```MainSlurm.jl``` file. The following breaks down the basic structure:
```Project.toml```: Contains the packages needed for the project. Each entry specifies the package name and its version requirement
```Initialize.jl```: Contains `Pkg.instantiate()` which reads the `Project.toml` file and automatically creates a project environment with the specified dependencies and versions. `Pkg.instantiate()` will trigger package installation when run for the first time
```MainSlurm.jl```: Main file to run the simulation
```src```: Contains the source code for the project, including all relevant functions that are called upon in ```MainSlurm.jl```
```examples```: Contains the example quadratic field configuation, example input image and a bash file `gen_kspace.sh` that calls upon `MainSlurm.jl`
```RunTests.jl```: Compiles the code and runs tests on a gradient echo and spin echo sequence using files from ```tests```
```tests```: Contains test scripts called upon in `RunTests.jl`
```FinalImageSim```: Final file to be run once `MainSlurm.jl` has finished running, to create the image from the simulated k-space


## Hardware Requirements
This is a very computationally heavy code and has been tested on High Performance Computing Systems. Julia uses [`threads`](https://docs.julialang.org/en/v1/manual/multi-threading/) to speed up the calculation by splitting it on multiple cores. To run the example gen_kspace.sh file, you will need:
- A computer with at least 256GB of RAM
- A multi-core CPU

## License
1<sup>ST</sup> ORDER MRI SIMS is MIT licensed, as found in [LICENSE file](https://github.com/radhikatibrewala/1stOrderMRISims/LICENSE)

## Cite
If you use the 1<sup>ST</sup> ORDER MRI SIMS code in your research, please use the following BibTeX entry.

```

```