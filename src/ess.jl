###############################################################################
# effective sample size calculations
###############################################################################

# univariate case
function ess(x::AbstractVector{<:Real};kwargs...)
    λ  = var(x)                  # estimate variance under stationary distribution 
    σ² = mcvar(x;kwargs...)[1,1] # estimate asymptotic variance of MCMC estimate 
    return length(x)*λ/σ²
end
