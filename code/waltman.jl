using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
DDEProblemLibrary.importddeproblems()

using PGFPlotsX

using DelimitedFiles
using LinearAlgebra

function plot_waltman(alg; plotdensity = 10_000)
    # compute solution
    sol = solve(DDEProblemLibrary.prob_dde_RADAR5_waltman_5, MethodOfSteps(alg);
                reltol = 1e-6, abstol = [1e-21, 1e-21, 1e-21, 1e-21, 1e-9, 1e-9])
    ts = range(sol.t[1]; stop = sol.t[end], length = plotdensity)

    # load RADAR5 solution
    radar5 = readdlm(joinpath(@__DIR__, "..", "data", "waltman_radar5.out"))
    radar5_ts = radar5[:, 1]
    radar5_us = radar5[:, 2:7]

    # create log plot of first four components
    @pgf logaxis = SemiLogYAxis()
    for i in 1:4
        @pgf push!(logaxis, PlotInc({ thick }, Table(radar5_ts, radar5_us[:, i])),
                   PlotInc({ thick, dashed }, Table(ts, sol(ts; idxs = i).u)))
    end

    # create standard plot of last two components
    @pgf axis = Axis(
        PlotInc({ thick }, Table(radar5_ts, radar5_us[:, 5])),
        PlotInc({ thick, dashed }, Table(ts, sol(ts; idxs = 5).u)),
        PlotInc({ thick }, Table(radar5_ts, radar5_us[:, 6])),
        PlotInc({ thick, dashed }, Table(ts, sol(ts; idxs = 6).u)),
    )

    # create empty group plot
    @pgf gp = GroupPlot(
        {
            group_style = {
                group_size = "2 by 1",
                xlabels_at = "edge bottom",
                ylabels_at = "edge left",
                horizontal_sep = raw"0.1\textwidth",
            },
            xlabel = raw"$t$", xmin = ts[1], xmax = ts[end],
            width = raw"0.5\textwidth", height = raw"0.7\textheight",
            grid = "major",
            label_style = { font = raw"\footnotesize" },
            tick_label_style = { font = raw"\footnotesize" },
            legend_pos = "north east",
            no_marks
        },
        logaxis,
        axis
    )
end

# specify LinSolveFactorize due to https://github.com/YingboMa/RecursiveFactorization.jl/issues/4
pgfsave(joinpath(@__DIR__, "..", "figures", "generated", "waltman.tex"),
        plot_waltman(Rosenbrock23(linsolve = LinSolveFactorize(lu)));
        include_preamble = false)
