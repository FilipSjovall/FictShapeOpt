using LinearSolve, LinearSolvePardiso, SparseArrays, 
      StaticArrays,FerriteMeshParser, Ferrite, 
      IterativeSolvers, AlgebraicMultigrid, IncompleteLU    



function load_files()
    include("mesh_reader.jl")

    include("material.jl")

    include("element_routines.jl")

    include("fem.jl")

    include("assemElem.jl")

    include("assem.jl")

    include("sensitivities.jl")
end


init_hyper()

dh0 = deepcopy(dh);

using Profile, ProfileView
@profile begin
    g_hist,v_hist, OptIter = Optimize(dh)
end

Profile.view()

@profview Optimize(dh)