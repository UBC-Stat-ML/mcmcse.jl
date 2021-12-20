using mcmcse
using Test
using RCall

function Gibbs_sampler(mu1, mu2, a, b, rho, init, n)
    X = similar(init, n, 2)
    X[1, 1] = init[1]
    X[1, 2] = init[2] 
    for i in 2:n
        X[i, 1] = mu1 + (rho / b) * (X[i - 1, 2] - mu2) + sqrt(a - (rho^2) / b)*randn()
        X[i, 2] = mu2 + (rho / a) * (X[i, 1] - mu1) + sqrt(b - (rho^2) / a)*randn()
    end
    return X
end

include("test_batch_size.jl")
include("test_mcvar.jl")
include("test_ess.jl")
