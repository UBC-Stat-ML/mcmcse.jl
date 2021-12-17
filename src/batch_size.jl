###############################################################################
# automatic batch size determination
# code ported from the R package mcmcse:
#     James M. Flegal, John Hughes, Dootika Vats, Ning Dai, Kushagra Gupta,
#     and Uttiya Maji. (2021). mcmcse: Monte Carlo Standard Errors for
#     MCMC. R package version 1.5-0. Riverside, CA, and Kanpur, India.
###############################################################################

function batch_size(
    M::TM;
    method::String = "bm",
    fast::Bool     = true
) where {F <: AbstractFloat, TM <: AbstractMatrix{F}}
    n, p = size(M)
    @assert n > p
    omax = min(p, n - 1, floor(Int, 10 * log10(n)))
    rows = 1:n
    fast && (rows = last(rows, 50_000))
    xacf = autocov(M[rows, :], 0:omax)
    if any(xacf[1,:] .< eps(F))
        @warn "No variability observed in a component. Setting batch size to 1"
        b = 1
    else
        threshold = convert(F, 1.959963984540054 / sqrt(n))
        bf = min(batch_size(n, xacf, method, threshold), floor(Int, n / (p + 1)))
        n > 10 && (bf = min(bf, floor(Int, n/10)))
        b = floor(Int, bf)
    end
    return b
end

###############################################################################
# low level functions
##############################################################################

# Calculates batchsize from Gamma and Sigma from AR approximation. Uses the "optimal" batchsize parametric 
# method  from of Liu et al.
function batch_size(
    n::Int,
    xacf_mat::TM,
    method::String,
    threshold::F
) where {F <: AbstractFloat, TM <: AbstractMatrix{F}}

    max_order = size(xacf_mat, 1) - 1
    p         = size(xacf_mat, 2)
    num_sum   = zero(F)
    denom_sum = zero(F)
        
    for i in 1:p  
        Gamma, Sigma = arp_approx(xacf_mat[:,i], max_order, n, threshold)
        num_sum     += Gamma^2
        denom_sum   += Sigma^2
    end

    coeff   = (num_sum/denom_sum)^(1.0/3.0)
    b_const = n*(method == "bm" ? 1.0 : 1.5)
    b       = max(1, (b_const^(1.0/3.0)) * coeff)
    return b
end

# Calculate Gamma and Sigma from AR approximation
function arp_approx(
    xacf::TV,
    max_order::Int,
    n::Int,
    threshold::F
) where {F <: AbstractFloat, TV <: AbstractVector{F}}

    coefs, vars, order = ar_yw(max_order, xacf, xacf, threshold, n)
    spec  = vars / ((one(F) - sum(coefs))^2)
    foo   = zero(F)

    if order != 0
        for i in 1:order
            for k in 1:i
                foo = foo + coefs[i] * k * xacf[abs(k-i)+1] # convert k to floor?
            end
        end
        Gamma = 2.0*(foo + (spec-xacf[1])/2.0 * dot(1:order, coefs)) / (1.0 - sum(coefs))
    else
        Gamma = zero(F)
    end
    return (Gamma = Gamma, Sigma = spec)
end

# Wrapper function to calculate AR coefficients and variance from the output of Levinson algorithm (eureka).
function ar_yw(
    order_max::Int,
    r::TV,
    g::TV,
    threshold::F,
    n::Int
) where {F <: AbstractFloat, TV <: AbstractVector{F}}

    vars, coefs, order = eureka(order_max, r, g, threshold)

    if order > 0
        coef_vec = coefs[order, 1:order]
        var_pred = vars[order]
    else
        coef_vec = zeros(F, order)
        var_pred = r[1]
    end

    var_pred = var_pred * n / (n - (order+1))

    return (coefs = coef_vec, vars = var_pred, order = order)
end

# Using the levinson algorithm to calculate AR coefficients and variance.
# We have added a check which compares the coefficient to a threshold. If at any time a coefficient falls 
# below the threshold, we don't compute coefficients beyond that order (i.e. we don't always complete the 
# loop from 1 to order.max) and return. 
function eureka(
    order_max::Int,
    r::TV,
    g::TV,
    threshold::F
) where {F <: AbstractFloat, TV <: AbstractVector{F}}
    
    # allocate storage
    coefs = similar(r, (order_max, order_max))
    var   = similar(r, order_max)
    a     = similar(var)
    
    v = r[1]
    d = r[2]
    a[1] = 1.0
    coefs[1,1] = g[2]/v

    if abs(coefs[1,1]) <= threshold
        return (vars = var, coefs = coefs, order = 0)
    end

    q = coefs[1,1] * r[2]
    var[1] = (1.0 - coefs[1,1]*coefs[1,1]) * r[1]

    if order_max == 1
        return (vars = var, coefs = coefs, order = 1)
    end

    for l in 2:order_max
        a[l] = -d/v

        if l > 2 
            l1 = (l - 2.0) / 2.0
            l2 = l1 + 1.0

            for j in 2:l2
                hold = a[j]
                k = l-j+1
                a[j] = a[j] + a[l] * a[k]
                a[k] = a[k] + a[l] * hold
            end

            if (2.0*l1) != (l-2)
                a[l2 + 1] = a[l2 + 1] * (1.0 + a[l])
            end
        end

        v = v + a[l] * d
        coefs[l,l] = (g[l+1] - q)/v

        if abs(coefs[l,l]) <= threshold 
            return (vars = var, coefs = coefs, order = l-1)
        end

        # since coefficients will most likely be in decreasing order, coefs[l-1,l-1] will be the smallest.
        # Therefore checking the first occurence of where this is smaller than the threshold would be enough.

        for j in 1:(l-1)
            coefs[l,j] = coefs[l-1,j] + coefs[l,l] * a[l-j+1]
        end

        var[l] = var[l-1] * (1.0 - coefs[l,l]*coefs[l,l])

        if l == order_max
            return (vars = var, coefs = coefs, order = order_max)
        end

        d = 0.0
        q = 0.0

        for i in 1:l
            k = l-i+2
            d = d + a[i] * r[k]
            q = q + coefs[l,i] * r[k]
        end

    end

    return (vars = var, coefs = coefs, order = order_max)
end
