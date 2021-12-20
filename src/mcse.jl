###############################################################################
# asymptotic variance estimation methods
# code ported from the R package mcmcse:
#     James M. Flegal, John Hughes, Dootika Vats, Ning Dai, Kushagra Gupta,
#     and Uttiya Maji. (2021). mcmcse: Monte Carlo Standard Errors for
#     MCMC. R package version 1.5-0. Riverside, CA, and Kanpur, India.
###############################################################################

function mcse(
    x::TM;
    method::String      = "bm",
    r::Real             = 3,
    b::Int              = -1,
    size_method::String = "optimal"
) where {F <: AbstractFloat, TM <: AbstractMatrix{F}}
    
    # parse method
    if method == "lug"
        method = "bm"
        r = 3
    elseif method != "bm"
        throw("Method $method unknown or not yet implemented")
    end

    # check r
    if r > 5
        @warn "We recommend using r ≤ 5. r ∈ {1,2,3} are standard"
    elseif r < 1
        throw(ArgumentError("r cannot be less than 1."))
    end

    # check dimensions
    c    = F(0.5)
    n, p = size(x)
    if n < (p+1)
        throw(ArgumentError("sample size is insufficient for a Markov chain of this dimension"))
    end
    
    # set batch size
    if b < 0 || b >= n || floor(n/size) <= 1
        if size_method == "optimal"
            b = batch_size(x, method = method)
        elseif size_method == "sqroot"
            b = floor(Int, sqrt(n))
        elseif size_method == "cuberoot"
            b = floor(Int, cbrt(n))
        else
            throw(ArgumentError("invalid b=$b and unknown size_method $size_method"))
        end
    end

    a = floor(Int, n/b) # number of batches
    b == 1 && r != 1 && (r = 1)
    if b < 2r
        @info "estimated batch size is low, lugsail not required. setting r = 1"
        r = 1
    end

    # estimate variance
    if b == 1
        sig_mat = cov(x)
    else
        sig_mat = bm(x, b)
        if r > 1
            sig_mat = (1/(1-c))*sig_mat - (c/(1-c))*bm(x, floor(Int, b/r))
        end
    end

    # check if not posdef
    if !isposdef(sig_mat)
        @warn "Estimated matrix not positive definite. The chain might be highly correlated or very high dimensional. Consider increasing the sample size."
    end

    return sig_mat
end

# handle vectors as 1-column matrices
function mcse(x::TV;kwargs...) where {F <: AbstractFloat, TV <: AbstractVector{F}}
    mcse(reshape(x,(length(x), 1));kwargs...)
end

# convert non-float reals to Float64
function mcse(x::AbstractVecOrMat{<:Real};kwargs...)
    mcse(convert(Array{Float64},x);kwargs...)
end

###############################################################################
# low level functions
###############################################################################

# batch means
function bm(M::TM, b::Int) where {F <: AbstractFloat, TM <: AbstractMatrix{F}}
    n = size(M, 1)
    a = floor(Int, n/b)    # a batches of size b
    y = sum(M, dims=1)/n   # compute means using all the data
    
    # compute centered batch means
    cbms = [sum(M[(1 + b*i):(b*(i + 1)), :], dims=1)/b .- y for i in 0:(a-1)]
    B    = reduce(vcat, cbms) # convert to matrix
    (b/(a-1)) * (B' * B)   # return estimate of variance
end

