# test
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

using RCall

for rep in 1:100
    mvg = Gibbs_sampler(2, 50, 1, 1, 0.5, [2.0, 50.0], 10_000)
    @rput mvg
    @assert all([
        rcopy(R"mcmcse:::batchSize(mvg, fast = FALSE)") == batch_size(mvg, fast=false)
        rcopy(R"mcmcse:::batchSize(mvg, method = 'obm', fast = FALSE)") == batch_size(mvg, method = "obm")
        rcopy(R"mcmcse:::batchSize(mvg, method = 'bartlett', fast = FALSE)") == batch_size(mvg, method = "bartlett")
        rcopy(R"mcmcse:::batchSize(mvg, method = 'tukey', fast = FALSE)") == batch_size(mvg, method = "tukey")
    ])
end
@code_warntype batch_size(mvg)

# M=mvg
# fast=false
# F = Float64
# xacf_mat = xacf
# r = g = xacf = xacf_mat[:,2]
# order_max = max_order = size(xacf_mat,1)-1