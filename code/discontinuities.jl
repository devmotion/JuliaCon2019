using DiffEqProblemLibrary.DDEProblemLibrary
DDEProblemLibrary.importddeproblems()

using PGFPlotsX

function plot_2delays_analytic(ts)
    prob = DDEProblemLibrary.prob_dde_constant_2delays_scalar
    xmin, xmax = extrema(ts)

    # plot analytic solution
    @pgf ax = Axis(
        {
            xlabel = raw"$t$", ylabel = raw"$x(t)$",
            no_marks, grid = "major",
            xmin = xmin, xmax = xmax
        },
    PlotInc({ thick }, Table(; x = ts, y = prob.f.analytic.(prob.u0, prob.h, prob.p, ts)))
    )

    ax
end

pgfsave(joinpath(@__DIR__, "..", "figures", "generated", "discontinuities_2delays.tex"),
        plot_2delays_analytic(range(0; stop = 1, length = 10_000));
        include_preamble = false)

