# mcmcse.jl

This is a pure-Julia adaptation of the R package `mcmcse`, which you can find on [CRAN](https://cran.r-project.org/package=mcmcse) and on [Github](https://github.com/dvats/mcmcse).

## Installation

From the Julia REPL, enter package mode by typing `]`, then write
```julia
add https://github.com/UBC-Stat-ML/mcmcse.jl.git
``` 

## Details

Currently `mcmcse.jl` contains only the core funcionalities in `mcmcse`

- `mcvar` implements some of the `mcse`-like functions in `mcmcse` via multiple-dispatch 
- `batch_size` implements `mcmcse::batchSize()` for automatic batch size determination
- `ess` implements `mcmcse::ess()` for estimating the effective sample size 


## References

James M. Flegal, John Hughes, Dootika Vats, Ning Dai, Kushagra Gupta, and Uttiya Maji. (2021). mcmcse: Monte Carlo Standard Errors for MCMC. R package version 1.5-0. Riverside, CA, and Kanpur, India.


