# **Documentation** for `eds-classification.jl`

## Load dependencies

````julia
using CSV
using DataFrames
using StatsBase
using Plots
````

## Add mineral classification algorithms to the path

````julia
include("donarummo_classification.jl");
include("kandler_classification.jl");
include("kutuzov_classification.jl");
include("panta_classification.jl");
````

## Add visualization function to the path

````julia
include("categorical_histogram.jl")
````

## Import example EDS net intensity and atom percent data
The `donarummo_classification()` function requires EDS net intensity data, the `kandler_classification()` and `panta_classification()` functions require atom percent data, and the `kutuzov_classification()` function requires atom fraction data (although passing atom percent data works because the function automatically converts the data into fractional units).

````julia
net = CSV.read("eds_net_intensities.csv", DataFrame);
atp = CSV.read("eds_atom_percents.csv", DataFrame);
````

## Apply the mineral classification algorithms to the appropriate data sets

````julia
don = donarummo_classification(net);
kan = kandler_classification(atp);
kut = kutuzov_classification(atp);
pan = panta_classification(atp);
````

## Visualize the data
### Donarummo algorithm

````julia
categorical_histogram(don[:,1],
  sortby=:descend,
  xrotation=45
)
````

### Kandler algorithm

````julia
categorical_histogram(kan[:,1],
  sortby=:descend,
  xrotation=45
)
````

### Kutuzov algorithm

````julia
categorical_histogram(kut[:,1],
  sortby=:descend,
  xrotation=45
)
````

### Panta algorithm

````julia
categorical_histogram(pan[:,1],
  sortby=:descend,
  xrotation=45
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

