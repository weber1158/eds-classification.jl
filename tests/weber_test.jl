## Import packages / functions
using CategoricalArrays, CSV, DataFrames, Serialization, DecisionTree, MLJ
include("../algorithms/weber_classification.jl");


## Load example data 
net = CSV.read("./data/eds_net_intensities.csv", DataFrame); 


## Print second row to visualize the structure of the data 
display(net[2,:])


## Apply Donarummo classification algorithm the full table 
web = weber_classification(net);


## Display result for an observation (i.e., a row) in the input data 
observation = 2
println("Weber ID #", observation, ": ", web.predicted_class[observation])
println("Probability score: ", maximum(web.probabilities[observation,:]), "%")
display(web.probabilities[observation,:]) # Display all probability scores
