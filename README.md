# InformMe

*An information-theoretic tool for analyzing DNA methylation sequencing data.*


| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GarrettJenkinson.github.io/InformMe.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GarrettJenkinson.github.io/InformMe.jl/dev) | [![Build Status](https://travis-ci.com/GarrettJenkinson/InformMe.jl.svg?branch=master)](https://travis-ci.com/GarrettJenkinson/InformMe.jl) [![Codecov](https://codecov.io/gh/GarrettJenkinson/InformMe.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/GarrettJenkinson/InformMe.jl) [![Coveralls](https://coveralls.io/repos/github/GarrettJenkinson/InformMe.jl/badge.svg?branch=master)](https://coveralls.io/github/GarrettJenkinson/InformMe.jl?branch=master) |


A julia port of the informME software package. informME is written
in MATLAB, bash, C++, and R. By porting the code to julia, we
hope to streamline the software, and remove the requirement
that the user has a MATLAB software license.

!!! warning "WARNING"
    THIS PACKAGE IS A WORK IN PROGRESS. FOR NOW, YOU SHOULD USE THE EXTENSIVELY TESTED MATLAB VERSION AVAILABLE HERE: [https://github.com/GarrettJenkinson/informME](https://github.com/GarrettJenkinson/informME)


## Installation

Samtools must be installed and on your system path.

The only unregistered julia dependency is QuadDIRECT which can be installed
by typing ']' without quotes in the julia repl to bring up the pkg
prompt and typing
```julia
pkg> add https://github.com/timholy/QuadDIRECT.jl.git
```

For now, informME.jl is an unregistered package, and should be installed
from the git repository at the pkg prompt:

```julia
pkg> add https://github.com/GarrettJenkinson/InformMe.jl
```
or if you have a local clone of the repository:

```julia
pkg> add /path/to/local/repo/InformMe.jl
```
the installation can be tested at the pkg prompt by running

```julia
pkg> test InformMe
```

## Documentation

- [**STABLE**](https://GarrettJenkinson.github.io/InformMe.jl/stable) &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**](https://GarrettJenkinson.github.io/InformMe.jl/dev) &mdash; *documentation of the in-development version.*

## License

All code except for the maxent algorithm (contained in /src/maxent.jl)
is GPLv3 licensed. The maxent algorithm has the following licensing
information in its header:

```
# Permission was provided by the author, Ali Mahammad-Djafari,
# to modify and distribute this code with the informME method's
# code. The original MATLAB code has been ported by Garrett Jenkinson
# to the julia language. Contact the original author
# directly for use outside of informME.
#
# Author's website:
# http://djafari.free.fr/index.htm
#
# Author's paper with source code:
# Mohammad-Djafari A. (1992) A Matlab Program to Calculate the Maximum Entropy Distributions. In: Smith C.R., Erickson G.J., Neudorfer P.O. (eds) Maximum Entropy and Bayesian Methods. Fundamental Theories of Physics (An International Book Series on The Fundamental Theories of Physics: Their Clarification, Development and Application), vol 50. Springer, Dordrecht
```


## Citing

See `CITATION.bib` for the bibtex formatted citations. If using the package in an academic setting, please be sure to cite the relevant papers:

[1] Jenkinson, G., Pujadas, E., Goutsias, J., and Feinberg, A.P. (2017), Potential energy landscapes identify the information-theoretic nature of the epigenome, Nature Genetics, 49: 719-729.

[2] Jenkinson, G., Abante, J., Feinberg, A.P., and Goutsias, J. (2018), An information-theoretic approach to the modeling and analysis of whole-genome bisulfite sequencing data, BMC Bioinformatics, 19:87, https://doi.org/10.1186/s12859-018-2086-5.

[3] Jenkinson, G., Abante, J., Koldobskiy, M., Feinberg, A.P., and Goutsias, J. (2019), Ranking genomic features using an information-theoretic measure of epigenetic discordance, BMC Bioinformatics, 20:175, https://doi.org/10.1186/s12859-019-2777-6.
