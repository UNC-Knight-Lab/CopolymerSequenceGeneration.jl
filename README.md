# CopolymerSequenceGeneration

## Description

Hello! The CopolymerSequenceGeneration package can be used to generate sequences of polymer segments where monomers incorporate statistically into a chain. The following code uses a kinetic Monte Carlo strategy to simulate Mayo-Lewis dynamics of chain growth. Only the feed ratios for each chain or chain segment as well as theoretical or experimental reactivity ratios for each of the monomer classes are currently needed. The v.1.0.0 package supports a binary monomer system. 

## Installation and usage

A version of Julia >= 1.6 is required to run this package which can be downloaded [here](https://julialang.org/). The CopolymerSequenceGeneration package can be downloaded and installed for use with:

```pkg > add https://github.com/UNC-Knight-Lab/CopolymerSequenceGeneration.jl```

The package can be called at the start of a Julia code with:

```using CopolymerSequenceGeneration```

Feed ratios are provided to the simulation as an input CSV file. and the number of chains to be simulated are both given as inputs using the following function for simulation:

```run_seq(feed_ratio_matrix, number_of_chains)```

A sample code snippet and feed ratios for a triblock copolymer are provided for user reference. Data can also be exported automatically to the file path where the feed ratio CSV was found.
