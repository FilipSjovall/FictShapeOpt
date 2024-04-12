using Pkg
Pkg.activate()
using Mortar2D, ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Printf, JLD2, Statistics, AlgebraicMultigrid
using CairoMakie #, Plots
#
include("..//mesh_reader.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")

# - - - - - - - - - - - - #
# Plot objective function #
# - - - - - - - - - - - - #
# Load results
@load "results//seal//p = 1/p=1_packning.jld2"
f = Figure(size = (360,288), fontsize=  12)
ax = Axis(f[1, 1],
    xgridvisible = false,
    ygridvisible = false,
    title  = L"\text{Objective function}", # L enables LaTeX strings
    xlabel = L"\text{Iteration}", #
    ylabel = L"g_0", #
    xtickalign = 1,  # Ticks inwards
    ytickalign = 1   # # Ticks inwards
)
lines!(1:true_iteration,g_hist[1:true_iteration], color = :blue)
f
# - - - - - - - - - - - - - - - - - - #
# Extract and plot contact tractions  #
# - - - - - - - - - - - - - - - - - - #
function plotTraction()
    traction = ExtractContactTraction(a, ε, coord)
    X_c = []
    tract = []
    for (key, val) ∈ traction
        append!(X_c, coord[key, 1])
        append!(tract, val)
    end
    ϵᵢⱼₖ = sortperm(X_c)
    tract = tract[ϵᵢⱼₖ]
    X_c = X_c[ϵᵢⱼₖ]
    return X_c, tract
end
X_c,traction = plotTraction()
with_theme(theme_latexfonts()) do
f = Figure()
    ax = Axis(f[1, 1],
        xgridvisible = false,
        ygridvisible = false,
        title  = "Contact traction",
        xlabel = "X-coordinate [m]",
        ylabel = L"$\lambda$ [N/m^2]",
        # Align x-axis ticks inside the plot
        xtickalign = 1,
        # Align y-axis ticks inside the plot
        ytickalign = 1
    )
    #scatterlines!(Float64.(X_c), Float64.(traction), .-, color = :red, markercolor = :red)
    scatterlines!(Float64.(X_c), Float64.(traction), color = :red)
    f
end


# Heatmaps
using FerriteViz
using GLMakie
a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
plotter = FerriteViz.MakiePlotter(dh,a);
FerriteViz.solutionplot(plotter,colormap=:turbo)
