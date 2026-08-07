function categorical_histogram(data_col; sortby=:label, normalization=:count, kwargs...)
#
# For plotting categorical data in a histogram-like fashion
#
# Required dependencies: 
#  - StatsBase.jl 
#  - Plots.jl
#
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
