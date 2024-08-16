using Pkg
Pkg.activate()
using ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Printf, JLD2, Statistics, AlgebraicMultigrid
using CairoMakie #, Plots926
#
begin
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
end
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
h_cm  = 8
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 1200 # dpi

reso = (w_cm * px_per_cm / width)

#@load "results//seal//v5//packning.jld2"
@load "results//seal//v6(byt_till_denna)//packning.jld2"
begin
    f = Figure( resolution = (width,height), fontsize = 12,font="CMU", px_per_unit = reso)
    ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Objective function $|f|$ [N$^3$]",
           limits = (0, true_iteration, 4, -g_hist[true_iteration]*1.25),
           leftspinecolor = :blue,
           ylabelcolor = :blue,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
           xgridvisible = false, ygridvisible = false,
           ylabel = "Contact force constraint", xlabel = L"\text{Iteration}",
           limits = (0, true_iteration, -0.5, 0.5),
           rightspinecolor = :red,
           leftspinecolor = :blue,
           ylabelcolor = :red,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
    lines!(ax1,1:true_iteration,-g_hist[1:true_iteration], color = :blue )
    lines!(ax2,1:true_iteration,au_hist[1:true_iteration], color = :red )
    f
    Makie.save("optimization_history_seal.pdf",f)
end

begin
    f = Figure( resolution = (width,height),font="CMU", fontsize = 12, px_per_unit = reso)
    ax = Axis(f[1, 1],
        xgridvisible = false,
        ygridvisible = false,
        xlabel = L"\text{Iteration}", #
        ylabelrotation = 0,
        ylabel = L"$f$ [N$^3$]", #
        xtickalign = 1,  # Ticks inwards
        ytickalign = 1,   # # Ticks inwards
        xticks = 0:50:200,
        topspinevisible = false,
        rightspinevisible = false,
        xminorticksvisible = false, yminorticksvisible = false,
        limits = (0, 120, 40.0, 85.0),
    )
    #@load "results//seal//v2//p = 1/packning.jld2" g_hist true_iteration
    #lines!(1:true_iteration,(abs.(g_hist[1:true_iteration])), color = :blue, label = L"p = 1")
    #@load "results//seal//v2//p = 3/packning.jld2" g_hist true_iteration
    @load "results//seal//v6(byt_till_denna)//packning.jld2" g_hist true_iteration
    lines!(1:true_iteration,(g_hist[1:true_iteration]), color = :green, label = L"p = 3")
    #f[1,2] = Legend(f, ax, framevisible = false, orientation = :vertical, tellwidth = true, tellheight = false)
    axislegend(ax, position = :rc, framevisible = false)
    f
    Makie.save("optimization_history_seal.pdf",f)
end

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #



cm_convert = 28.3465
w_cm  = 13# 8
h_cm  = 8# 13
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 1200 # dpi
reso = (w_cm * px_per_cm / width)

f = Figure( resolution = (width,height), fontsize = 12, px_per_unit = reso)
ax = Axis(f[1, 1],
    xgridvisible = false,
    ygridvisible = false,
    #title  = L"\text{Contact traction λ}", # L enables LaTeX strings
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

begin
    # parameters that should have been saved?
    Δ          = -0.025
    ε          = 1e5
    nloadsteps = 10
    mp₁        = [180 80].*1e3     # [K G]
    mp₂        = [2.5 0.1].*1e3    #
    # # # # # # # # #
    #@load "results//seal//v5//LabOpt.jld2" dh
    @load "results//seal//v6(byt_till_denna)//packning.jld2" dh
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
    σx, σy,τ,σᵛᵐ = StressExtract(dh, a, mp₁, mp₂)
    # Ny shit
    τ_c = ExtractContactTractionVec(a, ε, coord)
    traction = zeros(size(a))
    for (key,val) in τ_c
        dofs = register[key,:]
        traction[dofs] = val
    end
    #
    vtk_grid("results/seal/v5/p = 3-contact", dh) do vtkfile
        vtk_point_data(vtkfile, dh, a) # displacement field
        vtk_point_data(vtkfile, σx, "σx")
        vtk_point_data(vtkfile, σy, "σy")
        vtk_point_data(vtkfile, τ, "τ")
        vtk_point_data(vtkfile, σᵛᵐ, "σᵛᵐ")
        vtk_point_data(vtkfile, dh, traction, "traction") # ny shit
    end
    xc3,tract3 = plotTraction()
end
#begin
    #lines!(convert(Vector{Float32},xc1),convert(Vector{Float32},tract1), color = :blue,  label = L"p = 1", linestyle = :solid)
    #lines!(convert(Vector{Float32},xc2),convert(Vector{Float32},tract2), color = :red,   label = L"p = 2")
    lines!(convert(Vector{Float32},xc3),convert(Vector{Float32},tract3), color = :green, label = L"p = 3", linestyle = :solid)
    @load "initiellt_tryck.jld2" iX itract
    lines!(convert(Vector{Float32},iX),convert(Vector{Float32},itract), color = :black, label = L"\text{Initial}", linestyle = :dash)
    axislegend(ax, position = :rc, framevisible = false)
    #f[2, 1] = Legend(f, ax, orientation=:horizontal,framevisible = false,tellwidth = false, tellheight = true, nbanks=2, colgap=5, halign= :right, valign = :top)
    f
#end
Makie.save("traction_seal.pdf",f)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - #
# LSQ Pressure profile  #
# - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
cm_convert = 28.3465
w_cm  = 13# 8
h_cm  = 8# 13
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 1200 # dpi
reso = (w_cm * px_per_cm / width)

f = Figure( resolution = (width,height), fontsize = 12, px_per_unit = reso)
ax = Axis(f[1, 1],
    xgridvisible = false,
    ygridvisible = false,
    #title  = L"\text{Contact traction λ}", # L enables LaTeX strings
    xlabel = L"\text{Horizontal position [mm]}", #
    #ylabelrotation = 3π/2,
    ylabel = L"$λ$ [MPa]", #
    xtickalign = 1,  # Ticks inwards
    ytickalign = 1,   # # Ticks inwards
    topspinevisible = false,
    rightspinevisible = false,
    xminorticksvisible = true, yminorticksvisible = true,
    limits = (0.34, 0.52, 0.0, 100.0),
)
# run
begin
    Δ          = -0.025
    ε          = 1e5
    nloadsteps = 10
    mp₁        = [180 80].*1e3     # [K G]
    mp₂        = [2.5 0.1].*1e3    #
    @load "results//seal//lsq_seal_v1//packning.jld2" dh
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
    σx, σy,τ,σᵛᵐ = StressExtract(dh, a, mp₁, mp₂)
    # Ny shit
    τ_c = ExtractContactTractionVec(a, ε, coord)
    traction = zeros(size(a))
    for (key,val) in τ_c
        dofs = register[key,:]
        traction[dofs] = val
    end
    #
    xc3,tract3 = plotTraction()
    λ_target = ones(length(nₛ),1)
    for (i,node) in enumerate(nₛ)
        x = dh.grid.nodes[node].x[1]
        pmax = 50
        mid  = 0.5
        P    = 6
        width= 0.12
        λ_target[i] = pmax*exp( -( ((x-mid)^2) / width^2 )^P )
        #λ_target[i] = pmax*(1-3000*(x-mid)^4)# h(x)
    end
end
    scatter!(convert(Vector{Float64},xc3), vec(sort(λ_target,dims=1)), color= :red, marker = 'x', label = "Target") # Testa marker :xcross
    lines!(convert(Vector{Float64},xc3),convert(Vector{Float64},tract3), color = :green, label = "Optimized", linestyle = :solid)
    @load "results//seal//lsq_seal_v1//initiellt_tryck.jld2" iX itract
    lines!(convert(Vector{Float64},iX),convert(Vector{Float64},itract), color = :blue, label = "Initial", linestyle = :dash)
    axislegend(ax, position = :rt, framevisible = false, patchsize=(50,10))
    f
    Makie.save("LSQ_profile.pdf",f)
#
@load "results//seal//lsq_seal_v1//packning.jld2"
begin
    f = Figure( resolution = (width,height), fontsize = 12,font="CMU", px_per_unit = reso)
    ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Objective function $f$ [N/mm$^2$]",
           limits = (0, true_iteration, 0, 125),
           leftspinecolor = :blue,
           ylabelcolor = :blue,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
           xgridvisible = false, ygridvisible = false,
           ylabel = "Volume constraint", xlabel = L"\text{Iteration}",
           limits = (0, true_iteration, -0.015, 0.015),
           rightspinecolor = :red,
           leftspinecolor = :blue,
           ylabelcolor = :red,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
    lines!(ax1,1:true_iteration,g_hist[1:true_iteration], color = :blue )
    lines!(ax2,1:true_iteration,v_hist[1:true_iteration], color = :red )
    f
    Makie.save("optimization_history_seal_ex2.pdf",f)
end
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot optimization history: f and g₁ using separate y-axis #
# / specifically for Cylinder - Block problem               #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
using JLD2
@load "results//lunarc//cyl_konvergerad_ordentligt//OptimizationVariablesy.jld2"
using CairoMakie
set_theme!(theme_latexfonts())
cm_convert = 28.3465
w_cm  = 13
h_cm  = 10
width = w_cm*cm_convert
height= h_cm*cm_convert
px_per_cm = 600 # dpi
reso = w_cm * px_per_cm / width
f = Figure( resolution = (width,height), fontsize = 12,font="CMU", px_per_unit = reso)
ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Objective function $|f|$ [N]",
           limits = (0, 419, 4, -g_hist[419]*1.25),
           leftspinecolor = :blue,
           ylabelcolor = :blue,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Volume constraint $g_1$", xlabel = L"\text{Iteration}",
           limits = (0, 419, v_hist[1], 0.5),
           rightspinecolor = :red,
           leftspinecolor = :blue,
           ylabelcolor = :red,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
lines!(ax1,1:419,-g_hist[1:419], color = :blue )
lines!(ax2,1:419,v_hist[1:419], color = :red )
f
Makie.save("optimization_history.pdf",f)

## plot traction in paraview
begin
    @load "results//lunarc//flat | volume0.9//OptimizationVariablesy2_filter.jld2"
    coord, enod = getTopology(dh);
    n_bot = getnodeset(dh.grid,"n_bot")
    n_top = getnodeset(dh.grid,"n_top")
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
    register = getNodeDofs(dh)
    #τ_c = ExtractContactTractionVec(a,ε,coord)
    τ_c = ExtractContactTractionVec(Ψ.*0, μ, coord) # coord räknas från dh som har x = X₀ + Ψ, vi vill alltså inte lägga på Ψ igen!
    traction = zeros(size(a))
    for (key,val) in τ_c
        dofs = register[key,:]
        traction[dofs] = val
    end
    vtk_grid("results/cylinder_traction_🚜", dh) do vtkfile
                vtk_point_data(vtkfile, dh, 0*Ψ) # displacement field
                vtk_point_data(vtkfile, dh, traction, "traction")
    end
end

## plot traction after 1 iteration paraview
begin
    #@load "results//lunarc//flat_first_design_updates//OptimizationVariablesy2.jld2"
    @load "results//lunarc//flat_first_design_updates//OptimizationVariablesy2_nofilter.jld2"
    coord, enod = getTopology(dh);
    n_bot = getnodeset(dh.grid,"n_bot")
    n_top = getnodeset(dh.grid,"n_top")
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
    register = getNodeDofs(dh)
    #τ_c = ExtractContactTractionVec(a,ε,coord)
    τ_c = ExtractContactTractionVec(Ψ, μ, coord)
    traction = zeros(size(a))
    for (key,val) in τ_c
        dofs = register[key,:]
        traction[dofs] = val
    end
    vtk_grid("results/cylinder_traction_🚜_iter1", dh0) do vtkfile
                vtk_point_data(vtkfile, dh0, Ψ) # displacement field
                vtk_point_data(vtkfile, dh0, traction, "traction")
    end
end


# # # # # # # # # # # # # # # # # # # # # # # # #
#  Plot heat maps using FerriteViz and GLMakie  #
# # # # # # # # # # # # # # # # # # # # # # # # #
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
