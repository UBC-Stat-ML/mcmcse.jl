@testset "batch_size" begin
    pass = true
    rep  = 0
    while pass && rep<100
        rep += 1
        mvg = Gibbs_sampler(2, 50, 1, 1, 0.5, [2.0, 50.0], 10_000)
        @rput mvg
        pass = all([
            rcopy(R"mcmcse:::batchSize(mvg, fast = FALSE)") == batch_size(mvg, fast=false)
            rcopy(R"mcmcse:::batchSize(mvg, method = 'obm', fast = FALSE)") == batch_size(mvg, method = "obm")
            rcopy(R"mcmcse:::batchSize(mvg, method = 'bartlett', fast = FALSE)") == batch_size(mvg, method = "bartlett")
            rcopy(R"mcmcse:::batchSize(mvg, method = 'tukey', fast = FALSE)") == batch_size(mvg, method = "tukey")
        ])
    end
    @test pass
end