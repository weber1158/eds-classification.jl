using CSV, DataFrames, DecisionTree, MLJ, CategoricalArrays, Statistics, Serialization 


## Load training data
data = CSV.read("./data/model_training_data_balanced.csv", DataFrame) 
data_no_missing = dropmissing(data) 


## Define labels
mineral_cat = CategoricalArray(data_no_missing.Mineral) 
class_levels = levels(mineral_cat) 
labels_int = Int.(mineral_cat.refs) 


## Extract training features
features = select(data_no_missing, Not(:Mineral)) 


## 15% holdout for validation
(Xremainder, Xvalidation), (yremainder, yvalidation) = 
    partition((features, labels_int), 0.85,
        rng=10, multi=true, stratify=labels_int) 


## Split remaining data into 70/30 training/testing sets
(Xtrain, Xtest), (ytrain, ytest) = 
    partition((Xremainder, yremainder), 0.7,
        rng=10, multi=true, stratify=yremainder)


## Convert to a DataFrame 
Xtrain_features = Matrix{Float64}(Xtrain) 
Xtest_features = Matrix{Float64}(Xtest) 
Xvalidation_features = Matrix{Float64}(Xvalidation)
ytrain_targets = ytrain 
ytest_targets = ytest 
yvalidation_targets = yvalidation 


## Train machine learning model 
model = build_forest(ytrain_targets, Xtrain_features)
train_predictions = apply_forest(model, Xtrain_features) 
test_predictions = apply_forest(model, Xtest_features)
validation_predictions = apply_forest(model, Xvalidation_features) 
 

## Evaluate accuracy 
train_accuracy = mean(train_predictions .== ytrain_targets) 
test_accuracy = mean(test_predictions .== ytest_targets) 
validation_accuracy = mean(validation_predictions .== yvalidation_targets) 

println("Train accuracy:      ", round(train_accuracy * 100, digits=3), "%")
println("Test accuracy:       ", round(test_accuracy * 100, digits=3), "%") 
println("Validation accuracy: ", round(validation_accuracy * 100, digits=3), "%") 


## Confusion Matrix 
println("Test Confusion Matrix:") 
display(DecisionTree.confusion_matrix(ytest_targets, test_predictions)) 

println("\nValidation Confusion Matrix:")
display(DecisionTree.confusion_matrix(yvalidation_targets, validation_predictions)) 


## Save the model as a function for future use
serialize("./models/weber_forest_model_v2.jls", Dict("model" => model, "class_levels" => class_levels)) 
