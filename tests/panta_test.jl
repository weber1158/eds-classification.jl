## Import packages / functions
using CategoricalArrays, CSV, DataFrames
include("../algorithms/panta_classification.jl");


## Load example data 
atp = CSV.read("./data/eds_atom_percents.csv", DataFrame); 


## Print the second row to visualize the structure of the data 
display(atp[2,:])


## Apply Donarummo classification algorithm the full table 
pan = panta_classification(atp);


## Display result for an observation (i.e., a row) in the input data 
observation = 2
println("Panta ID #", observation, ": ", pan[observation, 1])