module test_angular_coverage
    include("simulation.jl")
    using LsqFit
    using Plots

    angs = collect(0:10:360) .* pi ./ 180
    n_Tests = length(angs)
    n_lip = 125
    iter = 3
    skip = 10
    n_steps = 100
    dataset = _ensemble!(n_Tests, n_steps, angs, n_lip; skip=skip, filename = joinpath(@__DIR__, "..\\data\\angular_coverage_test.jld"), iter=iter)

    step_vals = collect(1.0:skip:n_steps)
    tot_lip = iter * n_lip
    f_T(t, p) = p[1] .* t .+ 1
    fit_ps = []
    init_p = [rand()]

    for test in 1:n_Tests
        av_cluster_size_at_step = similar(step_vals)
        for step_index in 1:length(step_vals)
            data_at_step = dataset[test].sim_data[step_index]
            av_cluster_size_at_step[step_index] = tot_lip/(tot_lip - sum(data_at_step.size) + data_at_step.clusters)
        end
        fit = curve_fit(f_T, step_vals, av_cluster_size_at_step, init_p)
        push!(fit_ps, fit.param)
        #init_p = fit.param
    end

    Plots.plot(Int.(round.(angs .* 180 ./ pi)), reduce(hcat, fit_ps)'[:, 1], xlabel = "Angular Coverage (degrees)", ylabel = "Growth Rate", legend = false)
    savefig(joinpath(@__DIR__, "..\\figures\\growth_rate_vs_angular_coverage.png"))
end