module mcmcse

using StatsBase, Statistics, LinearAlgebra

export batch_size
include("batch_size.jl")

export mcvar
include("mcvar.jl")

export ess
include("ess.jl")

end
