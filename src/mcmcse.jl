module mcmcse

using StatsBase, Statistics, LinearAlgebra

export batch_size
include("batch_size.jl")

export mcse
include("mcse.jl")

end
