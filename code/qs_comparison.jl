using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
DDEProblemLibrary.importddeproblems()

# load QS problem
using DDEProblemLibrary: prob_dde_qs

using PGFPlotsX

function plot_qs_comparison()
    # time points at which solution is computed and plotted
    ts = range(prob_dde_qs.tspan[1]; stop = prob_dde_qs.tspan[2], length = 10_000)

    # create empty group plot
    @pgf gp = GroupPlot(
        {
            group_style = {
                group_size = "2 by 2",
                xlabels_at = "edge bottom",
                ylabels_at = "edge left",
                horizontal_sep = "1cm",
                vertical_sep = "1.25cm"
            },
            xlabel = raw"$t$", ylabel = raw"$A(t)$",
            x_unit = raw"\hour", y_unit = raw"\mol\per\litre",
            xmin = prob_dde_qs.tspan[1], xmax = prob_dde_qs.tspan[2],
            width = "5.5cm", height = "3.5cm",
            grid = "major",
            label_style = { font = raw"\footnotesize" },
            tick_label_style = { font = raw"\footnotesize" },
            title_style = { font = raw"\footnotesize\bfseries" },
        }
    )

    # for different ODE solvers
    for alg in (BS3(), Tsit5(), Rosenbrock23(), Rodas4())
        # compute solution of third component
        sol = solve(prob_dde_qs, MethodOfSteps(alg); save_idxs = 3)

        # plot AHL dynamics
        @pgf push!(
            gp,
            {
                title = string(nameof(typeof(alg)))
            },
            PlotInc(
                {
                    no_markers
                },
                Table(; x = ts, y = sol(ts))),
            PlotInc(
                {
                    only_marks, mark = "o"
                },
                Table(; x = sol.t, y = sol.u)))
    end

    gp
end

pgfsave(joinpath(@__DIR__, "..", "figures", "generated", "qs_comparison.tex"),
        plot_qs_comparison();
        include_preamble = false)
