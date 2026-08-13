"""
    categorical_histogram(data::CategoricalArray; sortby=:label, normalization=:count, kwargs...)

Histogram chart for categorical data

This function requires the StatsBase.jl and Plots.jl packages; it uses Plots.bar() to create the visualization.

# REQUIRED ARGUMENTS
- `data` : Vector of categorical data 

# OPTIONAL ARGUMENTS
- `sortby` : Order of the histogram bars. Options include :alphabetical, :descend, :ascend. The default value is :alphabetical

- `normalization` : Units of the y-axis. Options include :cdf, :count, :cumcount, :percentage, :probability. The default value is :count 

- `kwargs` : Optional arguments used by the Plots.bar() function.

# RETURNS
- A visualization using Plots.bar() 

"""
function categorical_histogram(data_col; sortby=:label, normalization=:count, kwargs...)
	counts = StatsBase.countmap(data_col);
	if (sortby == :alphabetical) | (sortby == :alphabet)
		sorted_pairs = sort(collect(counts), by = x -> x[1]);
	elseif (sortby == :descend) | (sortby == :descending)
		sorted_pairs = sort(collect(counts), by = x -> -x[2]);
	elseif (sortby == :ascend) | (sortby == :ascending)
		sorted_pairs = sort(collect(counts), by = x -> x[2]);
	else
		error("Error in categorical_histogram() using the 'sortby=:label' argument. Must be set equal to :alphabetical, :descend, or :ascend. Omitting the 'sortby=:label' argument will sort the data in alphabetical order by default.");
	end

	labels = first.(sorted_pairs);
	vals = last.(sorted_pairs);

	raw_counts = [counts[l] for l in labels];
	n = sum(raw_counts);
	
	if (normalization == :count) | (normalization == :counts)
		vals = raw_counts;
	elseif (normalization == :probability)
		vals = raw_counts ./ n;
	elseif (normalization == :percentage)
		vals = raw_counts ./ n .* 100;
	elseif (normalization == :cumcount)
		vals = cumsum(raw_counts);
	elseif (normalization == :cdf)
		vals = cumsum(raw_counts) ./ n;
	else
		error("Error in categorical_histogram() using the 'normalization=:count' argument. Must be set equal to :cdf, :count, :cumcount, :percentage, or :probability.")
	end
	
	Plots.bar(labels, vals; kwargs...)
end
