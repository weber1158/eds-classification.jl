## Import packages / functions
using CSV, DataFrames
include("../algorithms/donarummo_classification.jl");


## Load example data 
net = CSV.read("./data/eds_net_intensities.csv", DataFrame); 


## Print the second row to visualize the structure of the data 
display(net[2,:])


## Apply Donarummo classification algorithm the full table 
don = donarummo_classification(net);


## Display result for an observation (i.e., a row) in the input data 
observation = 2
println("Donarummo ID #", observation, ": ", don[observation, 1])
