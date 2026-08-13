"""
    weber_classification(data::DataFrame)
 
Random forest model for mineral identification

Takes a DataFrame of EDS net intensities for the elements Na, Mg, Al, Si, P, K, Ca, Ti, and Fe, then converts the DataFrame into an array of 23 elemental ratios to match the training features of the random forest model. The random forest model was trained on EDS net intensity data from 18 reference minerals acquired from the Smithsonian Institution. 

# REQUIRED ARGUMENTS 
- `data` :: DataFrame containing a column for each of the following elements: Na, Mg, Al, Si, P, K, Ca, Ti, and Fe. The name of each column may be the full element name or its abbreviation. For instance, "Silicon" and "Si" are valid table variable names. Both the American and British spelling of "Aluminum" ("Aluminium") are also valid. Capitalization is not required but spelling is paramount. The values in the table should represent the measured EDS net intensity for each element.

# RETURNS 
A `NamedTuple`: (predicited_class, probabilities). If you execute `result = weber_classification(data)`, then you can access the mineralogy of each row in `data` by executing `result.predicited_class`. You can also access a `DataFrame` of probability scores based on the machine learning model by executing `result.probabilities`. This is an extremely useful feature for gauging uncertainty! None of the other mineral classification algorithms this feature.

"""
function weber_classification(input_data::DataFrame; 
    model_path=joinpath(@__DIR__, "..", "models", "weber_forest_model_v2.jls"))
    #
    ### Convert input table into an appropriate DataFrame 
    #
    # Extract columns
    Na = input_data[:,:Na]
    Mg = input_data[:,:Mg];
    Al = input_data[:,:Al];
    Si = input_data[:,:Si];
    P  = input_data[:,:P];
    K  = input_data[:,:K];
    Ca = input_data[:,:Ca];
    Ti = input_data[:,:Ti];
    Fe = input_data[:,:Fe];
    # Create new columns
    Esums = Na.+Mg.+Al.+Si.+P.+K.+Ca.+Ti.+Fe;
    c01 = (Mg.+Fe)./Al; c02 = (Mg.+Fe)./Si; c03 = (Ca.+Na)./Al; c04 = K./(K.+Na.+Ca);
    c05 = K./(Al.+Si);  c06 = Al./Si;       c07 = Fe./Si;       c08 = Ca./Si;
    c09 = K./Si;        c10 = K./Al;        c11 = Ca./Na;       c12 = P./Ca;               
    c13 = Ti./Fe;       c14 = Mg./Al;       c15 = Na./Esums;    c16 = Mg./Esums;    
    c17 = Al./Esums;    c18 = Si./Esums;    c19 = P./Esums;     c20 = K./Esums;     
    c21 = Ca./Esums;    c22 = Ti./Esums;    c23 = Fe./Esums;
    # Concatenate into an array
    data_array = [c01 c02 c03 c04 c05 c06 c07 c08 c09 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23];
    
    # Load model
    saved = deserialize(model_path)
    model, class_levels = saved["model"], saved["class_levels"]
    # Apply model to input data
    class_codes = 1:length(class_levels)  # matches how CategoricalArray.refs encodes classes
    predicted_code = apply_forest(model, data_array)
    probabilities  = apply_forest_proba(model, data_array, class_codes) .* 100
 
    return (
        predicted_class = class_levels[predicted_code],
        probabilities   = DataFrame(probabilities, Symbol.(class_levels))
    )
end