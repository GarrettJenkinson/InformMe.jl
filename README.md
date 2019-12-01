# InformMe

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GarrettJenkinson.github.io/InformMe.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GarrettJenkinson.github.io/InformMe.jl/dev)
[![Build Status](https://travis-ci.com/GarrettJenkinson/InformMe.jl.svg?branch=master)](https://travis-ci.com/GarrettJenkinson/InformMe.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/GarrettJenkinson/InformMe.jl?svg=true)](https://ci.appveyor.com/project/GarrettJenkinson/InformMe-jl)
[![Codecov](https://codecov.io/gh/GarrettJenkinson/InformMe.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/GarrettJenkinson/InformMe.jl)
[![Coveralls](https://coveralls.io/repos/github/GarrettJenkinson/InformMe.jl/badge.svg?branch=master)](https://coveralls.io/github/GarrettJenkinson/InformMe.jl?branch=master)


A julia port of the informME software package. informME is written
in MATLAB, bash, C++, and R. By porting the code to julia, we
hope to streamline the software, and remove the requirement
that the user has a MATLAB software license.

WARNING: THIS PACKAGE IS A WORK IN PROGRESS. FOR NOW, YOU SHOULD USE
THE EXTENSIVELY TESTED MATLAB VERSION AVAILABLE HERE:

https://github.com/GarrettJenkinson/informME


# Installation

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

# License
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

See `CITATION.bib` for the relevant reference(s).
