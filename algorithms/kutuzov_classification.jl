"""
    kutuzov_classification(data::DataFrame)

Mineral classification scheme from Kutuzov et al. (2026)

Takes a DataFrame of EDS atom percent data for the elements Na, Mg, Al, Si, Ca, Ti, and Fe, then it runs through a series of comparative criteria checks to assign a mineralogy to each row. More accurately, the Kutuzov algorithm deals with atom fractions, but this function has been written to normalize the elemental data to sum to 1 automatically. Note that the Kutuzov algorithm was originally created for the analysis of data collected with single particle ICP-TOFMS, not SEM-EDS. However, in general it also works well with EDS data.

# REQUIRED ARGUMENTS
- `data` :: DataFrame containing a column for each of the following elements: Na, Mg, Al, Si, Ca, Ti, and Fe. The name of each column may be the full element name or its abbreviation. For instance, "Silicon" and "Si" are valid table variable names. Both the American and British spelling of "Aluminum" ("Aluminium") are also valid. Capitalization is not required by spelling is paramount. The values in the table should represent the measured atom percent for each element.

# RETURNS
A `DataFrame` containing the mineral assignments for each row in the input DataFrame. The list of possible mineral assignments is given below:

1. Phyllosilicate (clay minerals and mica)
2. Augite (clinopyroxene, high Ca)
3. Diopside (clinopyroxene, low Fe)
4. Pigeonite (orthopyroxene, low Ca)
5. High-Fe Hornblende (amphibole)
6. High-Mg Hornblende (amphibole)
7. Chlorite (clay mineral)
8. Kaolinite (clay mineral)
9. Albite (Feldspar, Na endmember)
10. Anorthite (Feldspar, Ca endmember)
11. Hypersthene (orthopyroxene, a variety of enstatite; high Mg)
12. Quartz (SiO2)
13. Ca-dominant (calcite, apatite, fluorite, etc.)
14. Fe-dominant (hematite, goethite, magnetite, etc.)
15. Unknown (does not match any of the known criteria)

# REFERENCE
Kutuzov, S., Olesik, J., Lomax-Vogt, M., Carter, L., Lowry, G., Bland, G., Wielinski, J., Sullivan, R., & Gabrielli, P. (2026). Geochemical characterization of millions of individual atmospheric particles entrapped in Antarctic ice across the last glacial-interglacial transition. Scientific Reports, 16(1), 10556. https://doi.org/10.1038/s41598-026-45260-3

"""
function kutuzov_classification(df)
    # Check that input is a DataFrame
        @assert typeof(df) == DataFrame "Input must be a DataFrame"
    
        # Convert full-name table variables to abbreviations
	    # Na, Mg, Al, Si, Ca, Ti, and Fe
        element_names = ["Aluminum","Aluminium","Silicon","Iron","Sodium",
        "Magnesium","Calcium","Titanium"] # Vector{String}
        element_names = map(lowercase,element_names)
        varnames = names(df)
        varnames_lower = map(lowercase,varnames) # Vector{String}
        element_abbreviations = ["Al","Al","Si","Fe","Na","Mg","Ca","Ti"]
        for (n, abbreviation) in enumerate(element_abbreviations)
            if any(x -> contains(x, element_names[n]), varnames_lower)
                idx = contains.(varnames_lower, element_names[n])
                varnames[idx] .= abbreviation
            end
        end
        rename!(df,varnames) # Replace variable names in DataFrame
    
        # Ensure that all 7 of the necessary columns exist
        if sum(in.(varnames, Ref(["Na","Mg","Al","Si","Ca","Ti","Fe"]))) != 7
            error("Input must be a DataFrame containing a column for Na, Mg, Al, Si, Ca, Ti, and Fe. Check the spellings of the column names. Only full element names and abbreviations are valid.")
        end

		# Create a new DataFrame for just those elements
		df_elements = df[:, [:Na,:Mg,:Al,:Si,:Ca,:Ti,:Fe]]

		# Sum each row
		df_row_sums = sum.(eachrow(df_elements))

		# Overwrite df to be a DataFrame of atom fractions rather than atom percentages
		df = df_elements ./ df_row_sums

		# Define local functions
		function mineral_classification(DF)
			num_classifications = nrow(DF)
			classifications = fill("Unknown", num_classifications)
			  # Begin calling indexes
			# 01 Phyllosilicate
				idx = check_phyllosilicate(DF,classifications)
				classifications[idx] .= ["Phyllosilicate"]
			# 02 Augite
				idx = check_augite(DF,classifications)
				classifications[idx] .= ["Augite"]
			# 03 Diopside
				idx = check_diopside(DF,classifications)
				classifications[idx] .= ["Diopside"]
			# 04 Pigeonite
				idx = check_pigeonite(DF,classifications)
				classifications[idx] .= ["Pigeonite"]
			# 05 High Fe Hornblende
				idx = check_fe_hornblende(DF,classifications)
				classifications[idx] .= ["High-Fe Hornblende"]
			# 06 High Mg Hornblende
				idx = check_mg_hornblende(DF,classifications)
				classifications[idx] .= ["High-Mg Hornblende"]
			# 07 Chlorite
				idx = check_chlorite(DF,classifications)
				classifications[idx] .= ["Chlorite"]
			# 08 Kaolinite
				idx = check_kaolinite(DF,classifications)
				classifications[idx] .= ["Kaolinite"]
			# 09 Albite
				idx = check_albite(DF,classifications)
				classifications[idx] .= ["Albite"]
			# 10 Anorthite
				idx = check_anorthite(DF,classifications)
				classifications[idx] .= ["Anorthite"]
			# 11 Hypersthene
				idx = check_hypersthene(DF,classifications)
				classifications[idx] .= ["Hypersthene"]
			# 12 Quartz
				idx = check_quartz(DF,classifications)
				classifications[idx] .= ["Quartz"]
			# 13 Ca-dominant
				idx = check_ca_dominant(DF,classifications)
				classifications[idx] .= ["Ca-dominant"]
			# 14 Fe-dominant
				idx = check_fe_dominant(DF,classifications)
				classifications[idx] .= ["Fe-dominant"]

			return classifications
		end # end mineral classification function

		function check_phyllosilicate(T, classification_array)
			# Extract columns from the DataFrame
			Si = T[:, :Si]
			Al = T[:, :Al]
			Mg = T[:, :Mg]
			Fe = T[:, :Fe]
			# Logical conditions
			condition1 = (Si .> 0)
			condition2 = (Al .> 0)
			condition3 = (Mg .> 0) .| (Fe .> 0)
			condition4 = (((Al .+ Mg .+ Fe)./Si).>0.7) .& 
				(((Al .+ Mg .+ Fe)./Si).<4.0)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4 
			return index
		end
		function check_augite(T, classification_array)
			# Extract columns from the DataFrame
			Ca = T[:, :Ca]
			Mg = T[:, :Mg]
			Fe = T[:, :Fe]
			Si = T[:, :Si]
			# Logical conditions
			condition1 = (Ca.>0) .& (Mg.>0) .& (Fe.>0) .& (Si.>0)
			expression2 = (Ca.+Mg.+Fe)./Si
			condition2 = (expression2.>0.6) .& (expression2.<2.2)
			expression3 = (Mg.+Fe)./Si
			condition3 = (expression3.>0.275) .& (expression3.<1.6)
			expression4 = Ca./Si
			condition4 = (expression4.>0.3) .& (expression4.<2.0)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4
			return index
		end
		function check_diopside(T, classification_array)
			# Extract columns from the DataFrame
			Ca = T[:, :Ca]
			Mg = T[:, :Mg]
			Si = T[:, :Si]
			# Logical conditions
			condition1 = (Ca.>0) .& (Mg.>0) .& (Si.>0)
			expression2 = (Ca./Si)
			condition2 = (expression2.>0.25) .& (expression2.<1.0)
			expression3 = (Mg./Si)
			condition3 = (expression3.>0.25) .& (expression3.<1.0)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3
			return index
		end
		function check_pigeonite(T, classification_array)
			# Extract columns from the DataFrame
			Ca = T[:, :Ca]
			Mg = T[:, :Mg]
			Fe = T[:, :Fe]
			Si = T[:, :Si]
			# Logical conditions
			condition1 = (Ca.>0) .& (Mg.>0) .& (Fe.>0) .& (Si.>0)
			expression2 = (Ca./Si)
			condition2 = (expression2.>0.06) .& (expression2.<0.25)
			expression3 = (Mg./Si)
			condition3 = (expression3.>0.22) .& (expression3.<0.88)
			expression4 = (Fe./Si)
			condition4 = (expression4.>0.22) .& (expression4.<0.88)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4
			return index
		end
		function check_fe_hornblende(T, classification_array)
			# Extract columns from the DataFrame
			Ca = T[:, :Ca]
			Na = T[:, :Na]
			Mg = T[:, :Mg]
			Al = T[:, :Al]
			Fe = T[:, :Fe]
			Si = T[:, :Si]
			Ti = T[:, :Ti]
			# Logical conditions
			condition1 = (Ca.>0) .& (Al.>0) .& (Fe.>0) .& (Si.>0)
			expression2 = (Fe.+(2 .* Ti).+Al.+Mg)./Si
			condition2 = (expression2.>0.428) .& (expression2.<1.5)
			expression3 = ((2 .* Ca).+Na)./Si
			condition3 = (expression3.<1.0)
			condition4 = (Mg./Si).<0.143
			expression5 = (Al./Si)
			condition5 = (expression5.>0.143) .& (expression5.<0.571)
			condition6 = (Fe./Al).>Mg 
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4 .& 
				condition5 .& condition6
			return index
		end
		function check_mg_hornblende(T, classifiation_array)
			# Extract columns from the DataFrame
			Ca = T[:, :Ca]
			Mg = T[:, :Mg]
			Al = T[:, :Al]
			Fe = T[:, :Fe]
			Si = T[:, :Si]
			Ti = T[:, :Ti]
			# Logical conditions
			condition1 = (Ca.>0) .& (Al.>0) .& (Fe.>0) .& (Si.>0) .& (Mg.>0)
			expression2 = (Fe.+(2 .* Ti).+Al.+Mg)./Si
			condition2 = (expression2.>0.428) .& (expression2.<3.0)
			condition3 = (Fe./Si).<0.143
			expression4 = (Al./Si)
			condition4 = (expression4.>0.143) .& (expression4.<0.571)
			condition5 = (Mg.+Al).>Fe
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4 .& 
				condition5
			return index
		end
		function check_chlorite(T, classification_array)
			# Extract columns from the DataFrame
			Mg = T[:, :Mg]
			Al = T[:, :Al]
			Fe = T[:, :Fe]
			Si = T[:, :Si]
			# Logical conditions
			condition1 = (Al.>0) .& (Fe.>0) .& (Si.>0) .& (Mg.>0)
			expression2 = (Al.+Fe.+Mg)./Si
			condition2 = (expression2.>0.95) .& (expression2.<3.6)
			# Combine all conditions
			index = condition1 .& condition2
			return index
		end
		function check_kaolinite(T, classifiation_array)
			# Extract columns from the DataFrame
			Al = T[:, :Al]
			Si = T[:, :Si]
			# Logical conditions
			condition1 = (Al.>0) .& (Si.>0)
			expression2 = (Al./Si)
			condition2 = (expression2.>0.8) .& (expression2.<1.2)
			# Combine all conditions
			index = condition1 .& condition2
			return index
		end
		function check_albite(T, classification_array)
			# Extract columns from the DataFrame
			Al = T[:, :Al]
			Si = T[:, :Si]
			Na = T[:, :Na]
			Ca = T[:, :Ca]
			# Logical conditions
			condition1 = (Al.>0) .& (Si.>0) .& (Na.>0)
			expression2 = (Na./Si)
			condition2 = (expression2.>0.2) .& (expression2.<0.66)
			expression3 = (Al./Si)
			condition3 = (expression3.>0.15) .& (expression3.<1.33)
			expression4 = (Ca./Si)
			condition4 = (expression4.<0.15)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4
			return index
		end
		function check_anorthite(T, classifiation_array)
			# Extract columns from the DataFrame
			Al = T[:, :Al]
			Si = T[:, :Si]
			Na = T[:, :Na]
			Ca = T[:, :Ca]
			# Logical conditions
			condition1 = (Al.>0) .& (Si.>0) .& (Ca.>0)
			expression2 = (Ca./Si)
			condition2 = (expression2.>0.25) .& (expression2.<1.0)
			expression3 = (Al./Si)
			condition3 = (expression3.>0.5) .& (expression3.<2.0)
			expression4 = (Na./Si)
			condition4 = (expression4.<0.2)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4
			return index
		end
		function check_hypersthene(T, classification_array)
			# Extract columns from the DataFrame
			Na = T[:, :Na]
			Mg = T[:, :Mg]
			Al = T[:, :Al]
			Si = T[:, :Si]
			Ca = T[:, :Ca]
			Fe = T[:, :Fe]
			Ti = T[:, :Ti]
			# Logical conditions
			condition1 = (Si.>0)
			condition2 = (Mg.>0) .| (Fe.>0)
			expression3 = (Mg.+Fe)./Si
			condition3 = (expression3.>0.75) .& (expression3.<1.25)
			condition4 = (Ca.<0.1) .& (Al.<0.1) .& (Na.<0.1) .& 
				((Ti./Si).<0.1)
			# Combine all conditions
			index = condition1 .& condition2 .& condition3 .& condition4
			return index
		end
		function check_quartz(T, classification_array)
			condition1 = T[:, :Si] .> 0.9
			index = condition1
			return index
		end
		function check_ca_dominant(T, classification_array)
			condition1 = T[:, :Ca] .> 0.6
			index = condition1
			return index
		end
		function check_fe_dominant(T, classification_array)
			condition1 = T[:, :Fe] .> 0.6
			index = condition1
			return index
		end

	# Final step: Classify the mineralogy for each row in the input table
	minerals = DataFrame(Minerals = mineral_classification(df))
	return minerals
end