###############################################################################
# effective sample size calculations
###############################################################################

# univariate case
function ess(x::AbstractVector{<:Real};kwargs...)
    λ  = var(x)             # estimate variance under stationary distribution 
    σ² = mcvar(x;kwargs...) # estimate asymptotic variance of the MCMC estimate 
    return length(x)*λ/σ²
end
