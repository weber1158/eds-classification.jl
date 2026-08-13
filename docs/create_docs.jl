using Literate
Literate.markdown("./docs/eds_docs.jl", outputdir=pwd(), flavor = Literate.FranklinFlavor())
Literate.notebook("./docs/eds_docs.jl", outputdir=pwd(), flavor = Literate.FranklinFlavor())