using DelayDiffEq
using PGFPlotsX

# create Mackey-Glass equation
f_mackey_glass(u, h, p, t) = 2 * h(p, t - 2) / (1 + h(p, t - 2)^9.65) - u

"""
    plot_mackey_glass(u₀, ts)

Plot numerical solutions of the Mackey-Glass equation
```math
u'(t) = 2 u(t - 2) / (1 + u(t - 2)^9.65) - u(t)
```
for ``t ∈ ts`` with ``a ∈ as`` and ``u(t) = u₀`` for ``t ≤ min ts``.
"""
function plot_mackey_glass(u₀, ts)
    # obtain time points
    t₀, t₁ = extrema(ts)

    # define history function
    h(p, t) = u₀

    # define DDE problem
    prob = DDEProblem(f_mackey_glass, u₀, h, (t₀, t₁); constant_lags = [2])

    # define DDE solver
    alg = MethodOfSteps(Rosenbrock23())

    # compute solution
    sol = solve(prob, alg, saveat = ts)

    # plot solution
    @pgf ax = Axis(
        {
            xlabel = raw"$t$", ylabel = raw"$x(t)$",
            xmin = t₀, xmax = t₁,
            no_marks, grid = "major"
        },
        PlotInc(Table(; x = ts, y = sol.u))
    )

    ax
end

"""
    plot_mackey_glass_embedding(u₀, ts)

Plot time delay embeddings, i.e., plots of ``u(t)`` versus ``u(t - 2)``,
of the Mackey-Glass equation
```math
u'(t) = 2 u(t - 2) / (1 + u(t - 2)^9.65) - u(t)
```
for ``t ∈ ts`` with ``a ∈ as`` and ``u(t) = u₀`` for ``t ≤ min ts``.
"""
function plot_mackey_glass_embedding(u₀, ts)
    # obtain time points
    t₀, t₁ = extrema(ts)

    # compute time points for time delay embedding
    Ts = filter(x -> x ≥ t₀ + 2, ts)

    # define history function
    h(p, t) = u₀

    # define DDE problem
    prob = DDEProblem(f_mackey_glass, u₀, h, (t₀, t₁); constant_lags = [2])

    # define DDE solver
    alg = MethodOfSteps(Tsit5())

    # compute solution
    sol = solve(prob, alg)

    # plot solution
    @pgf ax = Axis(
        {
            xlabel = raw"$x(t)$", ylabel = raw"$x(t - \tau)$",
            no_marks, grid = "major"
        },
        PlotInc(Table(; x = sol(Ts), y = sol(Ts .- 2)))
    )

    ax
end

pgfsave(joinpath(@__DIR__, "..", "figures", "generated", "mackey_glass.tex"),
        plot_mackey_glass(0.5, range(0; stop = 600, length = 10_000));
        include_preamble = false)
pgfsave(joinpath(@__DIR__, "..", "figures", "generated", "mackey_glass_embedding.tex"),
        plot_mackey_glass_embedding(0.5, range(0; stop = 600, length = 10_000));
        include_preamble = false)
