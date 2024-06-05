using Pkg
Pkg.activate()
using ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Printf, JLD2, Statistics, AlgebraicMultigrid
using CairoMakie #, Plots
#
include("Contact//Mortar2D//Mortar2D.jl")
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
    xticks = 0:50:200,
    topspinevisible = false,
    rightspinevisible = false,
    xminorticksvisible = false, yminorticksvisible = false,
    limits = (0, 225, 1.0, 4.0),
)
@load "results//seal//v2//p = 1/packning.jld2" g_hist true_iteration
lines!(1:true_iteration,(g_hist[1:true_iteration]./g_hist[1]), color = :blue, label = L"p = 1")
@load "results//seal//v2//p = 2/packning.jld2" g_hist true_iteration
lines!(1:true_iteration,(g_hist[1:true_iteration]./g_hist[1]), color = :red, label = L"p = 2")
@load "results//seal//v2//p = 3/packning.jld2" g_hist true_iteration
lines!(1:true_iteration,(g_hist[1:true_iteration]./g_hist[1]), color = :green, label = L"p = 3")
#f[2, 1] = Legend(f, ax, L"\text{Exponent}", framevisible = true, orientation = :horizontal, tellwidth = false, tellheight = true)
f[2, 1] = Legend(f, ax, framevisible = false, orientation = :horizontal, tellwidth = false, tellheight = true)
f
Makie.save("optimization_history_seal.pdf",f)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# Extract and plot contact tractions  #
function plotTraction()
    traction = ExtractContactTraction(a, ε, coord)
    X_c   = []
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
cm_convert = 28.3465
w_cm  = 8
h_cm  = 13
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 1200 # dpi
reso = (w_cm * px_per_cm / width)

f = Figure( resolution = (width,height), fontsize = 12, px_per_unit = reso)
ax = Axis(f[1, 1],
    xgridvisible = false,
    ygridvisible = false,
    title  = L"\text{Contact traction λ}", # L enables LaTeX strings
    xlabel = L"\text{Horizontal position [mm]}", #
    #ylabelrotation = 3π/2,
    ylabel = L"$λ$ [MPa]", #
    xtickalign = 1,  # Ticks inwards
    ytickalign = 1,   # # Ticks inwards
    topspinevisible = false,
    rightspinevisible = false,
    xminorticksvisible = true, yminorticksvisible = true,
    limits = (0.35, 0.5, 1.0, 200.0),
)
#### params
Δ          = -0.025
ε          = 1e5
nloadsteps = 10
mp₁        = [180 80].*1e3     # [K G]
mp₂        = [2.5 0.1].*1e3    #
begin
    @load "results//seal//v2//p = 1/packning.jld2" dh
    coord, enod = getTopology(dh);
    n_bot = getnodeset(dh.grid,"n_bot")
    n_top = getnodeset(dh.grid,"n_top")
    n_sym = getnodeset(dh.grid,"n_sym")
    nₛ    = getnodeset(dh.grid,"nₛ")
    nₘ    = getnodeset(dh.grid,"nₘ")
    contact_dofs = getContactDofs(nₛ, nₘ)
    contact_nods = getContactNods(nₛ, nₘ)
    freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)
    Γs = getfaceset(dh.grid,"Γ_slave")
    Γm = getfaceset(dh.grid,"Γ_master")
    global order = Dict{Int64,Int64}()
    for (i, nod) ∈ enumerate(contact_nods)
        push!(order, nod => i)
    end
    t = 1
    a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
    xc1,tract1 = plotTraction()

    @load "results//seal//v2//p = 2/packning.jld2" dh
    coord, enod = getTopology(dh);
    n_bot = getnodeset(dh.grid,"n_bot")
    n_top = getnodeset(dh.grid,"n_top")
    n_sym = getnodeset(dh.grid,"n_sym")
    nₛ    = getnodeset(dh.grid,"nₛ")
    nₘ    = getnodeset(dh.grid,"nₘ")
    contact_dofs = getContactDofs(nₛ, nₘ)
    contact_nods = getContactNods(nₛ, nₘ)
    freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)
    Γs = getfaceset(dh.grid,"Γ_slave")
    Γm = getfaceset(dh.grid,"Γ_master")
    global order = Dict{Int64,Int64}()

    for (i, nod) ∈ enumerate(contact_nods)
        push!(order, nod => i)
    end
    t = 1
    a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
    xc2,tract2 = plotTraction()


    @load "results//seal//v2//p = 3/packning.jld2" dh
    coord, enod = getTopology(dh);
    n_bot = getnodeset(dh.grid,"n_bot")
    n_top = getnodeset(dh.grid,"n_top")
    n_sym = getnodeset(dh.grid,"n_sym")
    nₛ    = getnodeset(dh.grid,"nₛ")
    nₘ    = getnodeset(dh.grid,"nₘ")

    contact_dofs = getContactDofs(nₛ, nₘ)
    contact_nods = getContactNods(nₛ, nₘ)
    freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)
    Γs = getfaceset(dh.grid,"Γ_slave")
    Γm = getfaceset(dh.grid,"Γ_master")
    global order = Dict{Int64,Int64}()
    for (i, nod) ∈ enumerate(contact_nods)
        push!(order, nod => i)

    end
    t = 1
    a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
    xc3,tract3 = plotTraction()
   #return xc1, tract1, xc2, tract2, xc3, tract3
end

lines!(convert(Vector{Float32},xc1),convert(Vector{Float32},tract1), color = :blue,  label = L"p = 1")
lines!(convert(Vector{Float32},xc2),convert(Vector{Float32},tract2), color = :red,   label = L"p = 2")
lines!(convert(Vector{Float32},xc3),convert(Vector{Float32},tract3), color = :green, label = L"p = 3")
@load "initiellt_tryck" X_c tract
lines!(convert(Vector{Float32},X_c),convert(Vector{Float32},tract), color = :black, label = L"\text{Initial}")
f[2, 1] = Legend(f, ax, orientation=:horizontal,framevisible = false,tellwidth = false, tellheight = true, nbanks=2, colgap=5, halign= :right, valign = :top)
f
Makie.save("traction_seal.pdf",f)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot optimization history: f and g₁ using separate y-axis #
# / specifically for Cylinder - Block problem               #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
using JLD2
#@load "färdig_cyl.jld2"
#@load "results//Cylinder + Platta//Bäst cylinder + platta//färdig_cyl.jld2"
@load "results//lunarc//flat | volume0.9//OptimizationVariablesy.jld2"
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
           ylabel = L"Objective function $|f|$ [N]",
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




### Från Jakob
function default_theme2(; size = (140, 70), rasterize=600, font="CMU", fontsizemajor=8, fontsizeminor=6)
   rasterize !== false ? rasterize = rasterize / 72 : nothing
   theme = merge(
       Theme(
           fontsize=fontsizemajor,
           font=font,
           figure_padding=2,
           size = size .* 2.83465,
           CairoMakie=Attributes(
               px_per_unit=5.0,
               pt_per_unit=2.0
           ),
           Axis=Attributes(
               xtickalign=1,
               ytickalign=1,
               xticklabelsize=fontsizeminor,
               yticklabelsize=fontsizeminor,
               xgridcolor=RGBAf(0.96, 0.96, 0.96),
               ygridcolor=RGBAf(0.96, 0.96, 0.96)
           ),
           Legend=Attributes(
               labelsize=fontsizeminor,
               titlesize=fontsizemajor,
               titlefont=:regular,
               rowgap=0,
               titlegap=0
           ),
           Label=Attributes(
               fontsize=fontsizemajor,
               titlesize=fontsizemajor
           ),
           Scatter=Attributes(
               strokewidth = 0.4
           ),
           ScatterLines=Attributes(
               strokewidth = 0.4,
               cycle = Cycle([:color, :linestyle, :marker], covary = true)
           ),
           Heatmap=Attributes(
               rasterize=rasterize
           ),
           Colorbar = Attributes(
               labelsize=fontsizemajor,
               ticklabelsize=fontsizeminor
           ),
           Lines = Attributes(
               linewidth = 1,
               cycle = Cycle([:color, :linestyle], covary = true)
           )
       ),
       theme_latexfonts()
   )
   return theme
end
