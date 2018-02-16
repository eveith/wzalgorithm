# wzAlgorithm

## Introduction

This shared library contains optimization algorithms that were part of the
Winzent project. Currently, these are:

  - Standard Particle Swarm Optimization (SPSO) 2011
  - The Multipart Evolutionary Strategy (REvol)

REvol has been described in the following scientific publication:

Martin Rupper, Eric MSP Veith, and Bernd Steinbach. "An Evolutionary
Training Algorithm for Artificial Neural Networks with Dynamic Offspring
Spread and Implicit Gradient Information."

## Installation

wzAlgorithm uses CMake for building. It requires Boost.Random.

    mkdir build; cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE ..
    make
    make test && make install

## Usage

Both REvol and SPSO live in a class of their own. These classes provide a
`run` function that starts the actual algorithm. `run` takes an evaluator
functor that takes the current individual/particle and returns `true` if the
individual/particle satisfies the optimization criterion, or `false` if the
search should continue. 

Both algorithms have tuneable parameters. The header file contains the bulk of
the documentation in Doxygen format.

## License

GNU GPLv3.
