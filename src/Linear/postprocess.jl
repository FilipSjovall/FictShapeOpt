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

# Kombinera dessa i en bild..
#
#
using JLD2
using CairoMakie
set_theme!(theme_latexfonts())
cm_convert = 28.3465
w_cm  = 13
h_cm  = 13
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 1200 # dpi

reso = (w_cm * px_per_cm / width)

f = Figure( resolution = (width,height), fontsize = 12, px_per_unit = reso)
ax = Axis(f[1, 1],
    xgridvisible = false,
    ygridvisible = false,
    title  = L"\text{Optimization history}", # L enables LaTeX strings
    xlabel = L"\text{Iteration}", #
    #ylabelrotation = 3π/2,
    ylabel = L"$\frac{f}{f^0}$", #
    xtickalign = 1,  # Ticks inwards
    ytickalign = 1,   # # Ticks inwards
    topspinevisible = false,
    rightspinevisible = false,
    xminorticksvisible = true, yminorticksvisible = true,
    limits = (0, 110, 1.0, 2.5),
)
@load "results//seal//p = 1/p=1.jld2" g_hist true_iteration
lines!(1:true_iteration,(g_hist[1:true_iteration]./g_hist[1]), color = :blue, label = L"p = 1")
@load "results//seal//p = 2/p=2.jld2" g_hist true_iteration
lines!(1:true_iteration,(g_hist[1:true_iteration]./g_hist[1]), color = :red, label = L"p = 2")
@load "results//seal//p = 3/p=3.jld2" g_hist true_iteration
lines!(1:true_iteration,(g_hist[1:true_iteration]./g_hist[1]), color = :green, label = L"p = 3")
#f[2, 1] = Legend(f, ax, L"\text{Exponent}", framevisible = true, orientation = :horizontal, tellwidth = false, tellheight = true)
f[2, 1] = Legend(f, ax, framevisible = false, orientation = :horizontal, tellwidth = false, tellheight = true)
f
Makie.save("optimization_history_seal.pdf",f)
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot optimization history: f and g₁ using separate y-axis #
# / specifically for Cylinder - Block problem               #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
using JLD2
#@load "färdig_cyl.jld2"
@load "results//Cylinder + Platta//Bäst cylinder + platta//färdig_cyl.jld2"
using CairoMakie
set_theme!(theme_latexfonts())
cm_convert = 28.3465
w_cm  = 13
h_cm  = 10
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 600 # dpi
reso = w_cm * px_per_cm / width
#f = Figure(fontsize=12,pt_per_unit = 1, px_per_unit = reso, size = (3000,2400) )
f = Figure( resolution = (width,height), fontsize = 12, px_per_unit = reso)
ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Objective function $f$ [N]",
           limits = (0, 400, 4, -g_hist[true_iteration]*1.25),
           leftspinecolor = :blue,
           ylabelcolor = :blue,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Volume constraint $g_1$", xlabel = L"\text{Iteration}",
           limits = (0, 400, v_hist[1], 0.5),
           rightspinecolor = :red,
           leftspinecolor = :blue,
           ylabelcolor = :red,
           xminorticksvisible = true, yminorticksvisible = true,
           title = L"\text{Optimization history}",
           topspinevisible = false)
#ax3 = Axis(f[1,1], width=Relative(0.2), height=Relative(0.2), halign=0.1, valign=1.0, backgroundcolor=:white)
lines!(ax1,1:true_iteration,-g_hist[1:true_iteration], color = :blue )
lines!(ax2,1:true_iteration,v_hist[1:true_iteration], color = :red )
#lines!(ax3,1:20,-g_hist[1:20],color =:blue )
f
# Makie.save("optimization_history.svg",f)
Makie.save("optimization_history.pdf",f)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Plot heat maps using FerriteViz and GLMakie          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Heatmaps
# using FerriteViz
# using GLMakie
# a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
# plotter = FerriteViz.MakiePlotter(dh,a);
# FerriteViz.solutionplot(plotter,colormap=:jet)
