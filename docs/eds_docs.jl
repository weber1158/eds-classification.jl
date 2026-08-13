# # **Documentation** for `eds-classification.jl`

# ## Import pacakges and functions
using CSV, DataFrames, Serialization, DecisionTree, MLJ
using StatsBase, Plots
include("./algorithms/edsClassification.jl");
include("./algorithms/categorical_histogram.jl");

# ## Check docstrings for a function
# Replace the function name as needed.
println(@doc weber_classification)


# ## Load example EDS data sets
net = CSV.read("./data/eds_net_intensities.csv", DataFrame);
atp = CSV.read("./data/eds_atom_percents.csv", DataFrame);

# ## Apply algorithms to the appropriate data set
don = donarummo_classification(net);
kan = kandler_classification(atp);
kut = kutuzov_classification(atp);
pan = panta_classification(atp);
web = weber_classification(net);


# ## Visualize the mineral classification data
# ### Donarummo
categorical_histogram(don[:,1], sortby=:descend, xrotation=45)


# ### Kandler 
categorical_histogram(kan[:,1], sortby=:descend, xrotation=45)


# ### Kutuzov 
categorical_histogram(kut[:,1], sortby=:descend, xrotation=45)


# ### Panta 
categorical_histogram(pan[:,1], sortby=:descend, xrotation=45)


# ### Weber
categorical_histogram(web.predicted_class, sortby=:descend, xrotation=45)


# ## Printing outputs for the Weber algorithm
# ### Example 1 - Find the mineral ID for an observation (i.e., row) and its associated probability score
observation = 485;
println("Weber prediction for mineral #485: \n", web.predicted_class[observation], "\n\nProbabilities:\n")
display(web.probabilities[observation,:])

# ### Example 2 - Compare the outputs of several algorithms 
observation = 15;
println("Predictions for mineral #15:")
println("True mineral:    Albite")
println("Weber pred.:     ", web.predicted_class[observation])
println("Donarummo pred.: ", don[observation,1])
println("Kandler pred.    ", kan[observation,1])
println("Kutuzov pred.    ", kut[observation,1])
println("Panta pred.      ", pan[observation,1])
