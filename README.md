# Demonstration code: PDEOPT
Optimal control with L1-regularization, state and control box-constraints
This is a demonstration code for solving optimal control problems as described in the article "PDE-Constrained Optimization: Optimal control with L1-regularization, state and control box-constraints"

The necessary matricies and structures are are given in files, the code is written in Julia. To run the code one needs to install the following packages: AlgebraicMultigrid.jl  (https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl), IterativeSolvers.jl (https://github.com/JuliaMath/IterativeSolvers.jl), FileIO.jl (https://github.com/JuliaIO/FileIO.jl), JLD2.jl (https://github.com/JuliaIO/JLD2.jl). All are availible in Julia's package manager. 

To plot the results I highly recommend Makie (https://github.com/JuliaPlots/Makie.jl) - as I am unsure it works on all systems I give an alternative plotter which uses PyPlot (https://github.com/JuliaPy/PyPlot.jl) too.

## Usage guide (you know Julia)
Download the files to a folder and run main.jl. The problem parameters can be changed in the section control panel found after the function declarations in main. By default the code will merely calculate the result, if you wish to plot the results uncomment the corresponding using- and function call at the end of main.jl.

## Usage guide (you are new to Julia)
Download Julia (https://julialang.org/) and open the Julia console (i.e., run julia-1.0.x/bin/julia). You can then install the needed packaged by:
```julia
using Pkg
Pkg.add("IterativeSolvers")
Pkg.add("AlgebraicMultigrid")
Pkg.add("FileIO")
Pkg.add("JLD2")
```
Download the Julia code (main.jl and aux.jl) together with data-files and place them in the same folder. Navigate to said folder in the Julia console. Some basic navigation commands: 
```julia
pwd() # Display current path
>> "/home/usr"
cd("Documents") # Change directory 
pwd()
>> "/home/usr/Documents"
cd() # Return home
```
Once in the correct folder run the program by writing:
```julia
include("main.jl")
```
The main problem parameters can be changed in the section control panel below the function declarations in main.jl. By default the code will merely calculate the result, if you wish to plot the results uncomment the corresponding using and function call at the end of main.jl (after installing Makie or PyPlot). 
