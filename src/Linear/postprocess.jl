begin
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
    function axisarrows!(ax::Axis=current_axis(); kwargs...)
        points = [
                ax.xaxis.attributes[:endpoints].val[2],
                ax.yaxis.attributes[:endpoints].val[2],
            ]
        directions = [Vec2f(1, 0), Vec2f(0, 1)]
        arrows!(ax.parent.scene, points, directions; kwargs...)
    end
    using JLD2
    using CairoMakie
    set_theme!(theme_latexfonts())
    cm_convert = 28.3465
    w_cm  = 13
    h_cm  = 8
    width = w_cm*cm_convert
    height= h_cm*cm_convert
    px_per_cm = 1200 # dpi
    #
    reso = (w_cm * px_per_cm / width)
end
# - - - - - - - - - - - - #
# Plot objective function #
# - - - - - - - - - - - - #
# Load results

# Kombinera dessa i en bild..
#
#


# # # # # # # # # # # # # #
# Seal objective function #
# # # # # # # # # # # # # #
@load "results//seal//v6(byt_till_denna)//packning.jld2"
begin
    f = Figure( resolution = (width,height), fontsize = 12,font="CMU", px_per_unit = reso)
    ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Objective function $|f|$ [N$^3$]",
           limits = (0, true_iteration, 1/1e4, -g_hist[true_iteration]*1.25e-4),
           leftspinecolor = :blue,
           ylabelcolor = :blue,
           xminorticksvisible = true, yminorticksvisible = true,
           #ytickformat = values -> ["$(value)⋅10⁵" for value in values],
           ytickformat = "{:0.1f}×10⁴",
           topspinevisible = false)
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
           xgridvisible = false, ygridvisible = false,
           ylabel = "Volume constraint", xlabel = L"\text{Iteration}",
           limits = (0, true_iteration, -0.01, 0.01),
           rightspinecolor = :red,
           leftspinecolor = :blue,
           ylabelcolor = :red,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
    lines!(ax1,1:true_iteration,-g_hist[1:true_iteration]./1e4, color = :blue )
    lines!(ax2,1:true_iteration,v_hist[1:true_iteration], color = :red )
    f
    Makie.save("optimization_history_seal.pdf",f)
end

# # # # # # # # # # # # # # # # # # # # # # # #
# Seal objective function with no constraint  #
# # # # # # # # # # # # # # # # # # # # # # # #
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

begin
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
        xlabel = L"Horizontal position $x$ [mm]", #
        #ylabelrotation = 3π/2,
        ylabel = L"$λ$ [MPa]", #
        xtickalign = 1,  # Ticks inwards
        ytickalign = 1,   # # Ticks inwards
        topspinevisible = false,
        rightspinevisible = false,
        #xminorticksvisible = true, yminorticksvisible = true,
        limits = (0.35, 0.5, 1.0, 200.0),
    )
    @load "results//seal//v6(byt_till_denna)//LabOpt.jld2" X_c tract
    lines!(convert(Vector{Float32},X_c),convert(Vector{Float32},tract), color = :green, label = L"\text{Optimized}", linestyle = :solid)
    @load "results//seal//v6(byt_till_denna)//initiellt_tryck.jld2" iX itract
    lines!(convert(Vector{Float32},iX),convert(Vector{Float32},itract), color = :black, label = L"\text{Initial}", linestyle = :dash)
    axislegend(ax, position = :rt, framevisible = false)
    ax.yticks = [0,50,100,150]
    axisarrows!(arrowsize=10)
    hidespines!(ax, :t, :r)
    f
    Makie.save("traction_seal.pdf",f)
end
#end


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - #
# LSQ Pressure profile  #
# - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

# Plot LSQ - pressure profiles
begin
    function t_func(x)
        pmax = 60
        mid  = 0.5
        P    = 6
        width= 0.06
        return pmax*exp( -( ((x-mid)^2) / width^2 )^P )
    end
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
        xlabel = L"Horizontal position $x$ [mm]", #
        #ylabelrotation = 3π/2,
        ylabel = L"$λ$ [MPa]", #
        xtickalign = 1,  # Ticks inwards
        ytickalign = 1,   # # Ticks inwards
        topspinevisible = false,
        rightspinevisible = false,
        xminorticksvisible = false, yminorticksvisible = false,
        limits = (0.34, 0.52, 0.0, 100.0),
    )
    @load "results//seal//lsq_seal_v4//LabOpt.jld2" X_c tract λ_target
    lines!(convert(Vector{Float64},X_c),convert(Vector{Float64},tract), color = :green, label = "Optimized", linestyle = :solid)
    for (i,x) in enumerate(X_c)
        λ_target[i] = t_func(x)
    end
    scatter!(convert(Vector{Float64},X_c), vec(sort(λ_target,dims=1)), color= :red, marker = :xcross, markersize=10, label = "Target") # Testa marker :xcross
    @load "results//seal//lsq_seal_v4//initiellt_tryck.jld2" iX itract
    lines!(convert(Vector{Float64},iX),convert(Vector{Float64},itract), color = :blue, label = "Initial", linestyle = :dash)
    axislegend(ax, position = :rt, framevisible = false, patchsize=(50,10))
    # test
    ax.yticks = [0,25,50,75]
    axisarrows!(arrowsize=10)
    hidespines!(ax, :t, :r)
    # test
    f
    Makie.save("LSQ_profile.pdf",f)
end
# Objective function
@load "results//seal//lsq_seal_v4//LabOpt.jld2"
begin
    f = Figure( resolution = (width,height), fontsize = 12,font="CMU", px_per_unit = reso)
    ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
           xgridvisible = false, ygridvisible = false,
           ylabel = L"Objective function $f$ [N/mm$^2$]",
           limits = (0, true_iteration, 0, 20),
           leftspinecolor = :blue,
           ylabelcolor = :blue,
           xminorticksvisible = true, yminorticksvisible = true,
           topspinevisible = false)
    #ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
    #       xgridvisible = false, ygridvisible = false,
    #       ylabel = "Volume constraint", xlabel = L"\text{Iteration}",
    #       limits = (0, true_iteration, -0.015, 0.015),
    #       rightspinecolor = :red,
    #       leftspinecolor = :blue,
    #       ylabelcolor = :red,
    #       xminorticksvisible = true, yminorticksvisible = true,
    #       topspinevisible = false)
    lines!(ax1,1:true_iteration,g_hist[1:true_iteration], color = :blue )
    hidespines!(ax1, :t, :r)
    #lines!(ax2,1:true_iteration,v_hist[1:true_iteration], color = :red )
    f
    Makie.save("optimization_history_seal_ex2.pdf",f)
end
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot optimization history: f and g₁ using separate y-axis #
# / specifically for Cylinder - Block problem               #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
begin
    using JLD2
    @load "results//lunarc//samma_interferens_aktuell//OptimizationVariablesy.jld2"
    using CairoMakie
    set_theme!(theme_latexfonts())
    cm_convert = 28.3465
    w_cm  = 13
    h_cm  = 10
    width = w_cm*cm_convert
    height= h_cm*cm_convert
    px_per_cm = 600 # dpi
    n_iter = 417
    reso = w_cm * px_per_cm / width
    f = Figure( resolution = (width,height), fontsize = 12,font="CMU", px_per_unit = reso)
    ax1 = Axis(f[1, 1], yticklabelcolor = :blue,
            xgridvisible = false, ygridvisible = false,
            ylabel = L"Objective function $|f|$ [N]",
            limits = (1, n_iter, 2, 12),
            leftspinecolor = :blue,
            ylabelcolor = :blue,
            xminorticksvisible = true, yminorticksvisible = true,
            topspinevisible = false)
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right,
            xgridvisible = false, ygridvisible = false,
            ylabel = L"Volume constraint $g_1$", xlabel = L"\text{Iteration}",
            limits = (1, n_iter, v_hist[1], 0.5),
            rightspinecolor = :red,
            leftspinecolor = :blue,
            ylabelcolor = :red,
            xminorticksvisible = true, yminorticksvisible = true,
            topspinevisible = false)
    lines!(ax1,1:n_iter,-g_hist[1:n_iter]./2, color = :blue )
    lines!(ax2,1:n_iter,v_hist[1:n_iter], color = :red )
    f
    Makie.save("optimization_history.pdf",f)
end
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
    mp    = [175 80.769230769230759]
    Δ     = -0.05
    coord, enod = getTopology(dh);
    n_bot = getnodeset(dh.grid,"n_bot")
    n_top = getnodeset(dh.grid,"n_top")
    nₛ    = getnodeset(dh.grid,"nₛ")
    nₘ    = getnodeset(dh.grid,"nₘ")
    n_left  = getnodeset(dh.grid,"nₗ")
    n_right = getnodeset(dh.grid,"nᵣ")
    K       = create_sparsity_pattern(dh)
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

using Pkg
Pkg.activate()
using ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using Printf, JLD2, Statistics # AlgebraicMultigrid # SparseDiffToolss
using CairoMakie #, Plots926
set_theme!(theme_latexfonts())
# series = [5 25 50 75 100 200]
first_string = "results//seal//v7(finare)//tryck_at_"
last_string = ".jld2"
begin
    @load "results//seal//v7(finare)//initiellt_tryck_v2.jld2" iX itract
    f = plot_λ_series(iX, itract, 1)
    Makie.save("results//seal//v7(finare)//pressure_hist//initial_pressure_hist.pdf",f)
    for num in 5:5:410 #series
        load_str = first_string*string(num)*last_string
        @load load_str X_c tract
        f = plot_λ_series(X_c, tract,num)
        Makie.save("results//seal//v7(finare)//pressure_hist//pressure_hist_"*string(num)*".pdf",f)
    end
    @load "results//seal//v7(finare)//LabOpt_v2.jld2" X_c tract OptIter
    f = plot_λ_series(X_c, tract, OptIter+190)
    Makie.save("results//seal//v7(finare)//pressure_hist//pressure_hist_"*string(OptIter+190)*".pdf",f) # OptIter resettas fram till 200
    println(OptIter)
end
function plot_λ_series(x,y,num)
    function axisarrows!(ax::Axis=current_axis(); kwargs...)
        points = [
                ax.xaxis.attributes[:endpoints].val[2],
                ax.yaxis.attributes[:endpoints].val[2],
            ]
        directions = [Vec2f(1, 0), Vec2f(0, 1)]
        arrows!(ax.parent.scene, points, directions; kwargs...)
    end
    cm_convert = 28.3465
    w_cm  = 6 # 8
    h_cm  = 6 # 13
    width = w_cm*cm_convert
    height= h_cm*cm_convert
    px_per_cm = 1200 # dpi
    reso = (w_cm * px_per_cm / width)
    f = Figure( resolution = (width,height), fontsize = 8, px_per_unit = reso)
    ax = Axis(f[1, 1],
        xgridvisible = false,
        ygridvisible = false,
        title  = "Design iteration " * string(num), # L enables LaTeX strings
        xlabel = L"Horizontal position $x$ [mm]", #
        #ylabelrotation = 3π/2,
        ylabel = L"$λ$ [MPa]", #
        xtickalign = 1,  # Ticks inwards
        ytickalign = 1,   # # Ticks inwards
        topspinevisible = false,
        rightspinevisible = false,
        xminorticksvisible = false, yminorticksvisible = false,
        limits = (0.34, 0.52, -0.1, 220.0),
    )
    lines!(convert(Vector{Float64},x),convert(Vector{Float64},y), color = :blue, linestyle = :solid)
    #scatter!(x,y, color = :blue)
    # axislegend(ax, position = :rt, framevisible = false, patchsize=(50,10))
    # test
    # ax.yticks = [0,100,200]
    # ax.xticks = [0.35, 0.40, 0.50]
    axisarrows!(arrowsize=10)
    hidespines!(ax, :t, :r)
    f
    return f
end


#=
#
#
# # # # # # # # # # # # # # # # # # # # # # # # #
#  Plot heat maps using FerriteViz and GLMakie  #
# # # # # # # # # # # # # # # # # # # # # # # # #
# Heatmaps
# using FerriteViz
# using GLMakie
# a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
# plotter = FerriteViz.MakiePlotter(dh,a);
# FerriteViz.solutionplot(plotter,colormap=:jet)
#
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
#
#
=#
