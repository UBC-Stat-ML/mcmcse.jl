
using StatsBase
using LinearAlgebra

function mbm(M::AbstractMatrix, b::Int)
    n = size(M, 1)
    a = floor(Int, n/b)    # a batches of size b
    y = sum(M, dims=1)/n   # compute means using all the data
    # compute centered batch means
    cbms = [sum(M[(1 + b*i):(b*(i + 1)), :], dims=1)/b .- y for i in 0:(a-1)]
    B = reduce(vcat, cbms) # convert to matrix
    (b/(a-1)) * (B' * B)
end

using RCall

# test
M = randn(100,4)
b = ceil(Int, 30*rand())
R"""
varmat = mcmcse:::mbmC($M, $b)
"""
out_R = @rget varmat
out = mbm(M, b)
out ≈ out_R

