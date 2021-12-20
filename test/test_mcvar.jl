@testset "mcvar" begin
    pass = true
    rep  = 0
    while pass && rep<100
        rep += 1
        mvg  = Gibbs_sampler(2, 50, 1, 1, 0.5, [2.0, 50.0], 10_000)
        @rput mvg
        pass = rcopy(R"mcmcse:::mcse.multi(mvg)$cov") ≈ mcvar(mvg)
    end
    @test pass
end