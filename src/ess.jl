###############################################################################
# effective sample size calculations
# code ported from the R package mcmcse:
#     James M. Flegal, John Hughes, Dootika Vats, Ning Dai, Kushagra Gupta,
#     and Uttiya Maji. (2021). mcmcse: Monte Carlo Standard Errors for
#     MCMC. R package version 1.5-0. Riverside, CA, and Kanpur, India.
###############################################################################

# univariate case
function ess(x::AbstractVector{<:Real};kwargs...)
    lambda = var(x)                  # estimate variance under stationary distribution 
    σ²     = mcvar(x;kwargs...)[1,1] # estimate asymptotic variance of MCMC estimate 
    return length(x)*lambda/σ²
end
