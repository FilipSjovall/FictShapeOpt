using Zygote
using Base.Threads

function square(x)
    return x^2
end

function vectorSquare(x)
    vec = similar(x)
    @threads for (el, val) in enumerate(x)
        @inbounds vec[el] = square(x[el])
    end
    return vec
end

square(2)

vals = [1.0; 2.0; 3.0]

vectorSquare(vals)

Zygote.jacobian(x -> vectorSquare(x), vals)

gradient = ForwardDiff.gradient(vectorSquare, vals)

@threads for i in 1:10
    println("I'm thread $(threadid()), working on i=$i")
end

using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using SparseDiffTools
using Plots
using Printf
using JLD2
