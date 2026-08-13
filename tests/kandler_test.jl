## Import packages / functions
using CSV, DataFrames
include("../algorithms/kandler_classification.jl");


## Load example data 
atp = CSV.read("./data/eds_atom_percents.csv", DataFrame); 


## Print the second row to visualize the structure of the data 
display(atp[2,:])


## Apply Donarummo classification algorithm the full table 
kan = kandler_classification(atp);


## Display result for an observation (i.e., a row) in the input data 
observation = 2
println("Kandler ID #", observation, ": ", kan[observation, 1])
