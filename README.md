# Liposome Aggregation Simulation
## Introduction
This repository contains code for simulating the aggregation of liposomes in the presence of mRNA, using an Agent-Based Model (ABM). 

The advent of lipid-mediated mRNA delivery systems has dramatically changed the immunological landscape, particularly in showing great promise in personalized cancer immunotherapy. Given the importance and potential of this method of mRNA and drug delivery, it is imperative to understand the fundamental principles underlying these types of formulations for effective and precise medicine. 

Observations from real-time fluorescence imaging show that an assembly process governed by Smoluchowski kinetics underlies the formation of heterogeneous conglomerates spanning several orders of magnitude in size. The assembly process likely occurs in discrete steps after the mixing of mRNA and cationic liposomes in aqueous solution. Following the initial mRNA adsorption, the mRNA-coated liposomes self-assemble as they diffuse through their medium.

This simulation is written in Julia code and uses the `Agents.jl` package. The code models the movement of liposomes over time and observes the formation of aggregates through their interactions. Because the timescales of mRNA adsorption onto the liposomes is much faster than that of liposome aggregation, it is not modeled here. These simulation results can serve to verify our model of this process.

This work is based on our paper, *Multi-Step Assembly of an RNA-Liposome Nanoparticle Formulation Revealed by Real-Time, Single-Particle Quantitative Imaging*, published in *Advanced Science*. The full paper can be accessed at: [http://doi.org/10.1002/advs.202414305](http://doi.org/10.1002/advs.202414305).

## Code
The simulation code is organized into three main files:

- `model.jl`: Defines a model of the system, including liposomes and clusters. This includes the equations for liposome motion and interaction mechanics between liposomes.
- `simulation.jl`: Implements code to run the simulation once or as an ensemble and collect data from those simulations into custom data structures. It can also produce videos of the process (TODO).
- `analysis.jl`: Runs an experiment where the number of liposomes is kept constant and the angular coverage is varied. It saves the results of the simulation to `\data` and a figure of the growth rate vs angular coverage to `\figures` 

## Usage
  
To run the simulation, you'll first need to download the code. Here's how:  
1. Click the green "Code" button and select "Download ZIP." (alternatively you could clone this git repository if you are familiar with that)  
2. Extract the contents of the ZIP file to a location on your computer.  

Next, you'll need to run the simulation using Julia (If you don't have Julia installed, download and install it from [here](https://julialang.org/downloads/))  
1. Open a Julia terminal.  
2. Change the current directory to the location where you extracted the simulation code (e.g., by using the cd command).  
3. Run the following command:  
`include("src\\analysis.jl")`  

This will perform a single simulation and save the data object to `\data` and plots to `\figures`. 

## License
The code in this repository is released under the MIT License. Please see the LICENSE file for more information.
