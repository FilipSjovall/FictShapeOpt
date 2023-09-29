# Work in progress!
# Testa i OptimizationLinCircleTest
#
struct FEM
    K::SparseMatrixCSC{Float64,Int64}
    Kψ::SparseMatrixCSC{Float64,Int64}
    a::Vector{Float64}
    d::Vector{Float64}
    Ψ::Vector{Float64}
    Fᵢₙₜ::Vector{Float64}
    rc::Vector{Float64}
    Fₑₓₜ::Vector{Float64}
    Δa::Vector{Float64}
    res::Vector{Float64}
    bcdof::Vector{Int64}
    bcval::Vector{Float64}
    bcdof2::Vector{Int64}
    bcval2::Vector{Float64}
    mp::Matrix{Float64}
    mp₀::Matrix{Float64}
    t::Float64
end

struct Sensitivities
    g::Float64
    #global ∂rᵤ_∂x = similar(K) # behövs inte om vi har lokal funktion?
    #global dr_dd = similar(K) # behövs inte om vi har lokal funktion?
    #global ∂rψ_∂d = similar(K) # behövs inte om vi har lokal funktion?
    #global ∂g_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
    #global ∂g_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
    #global ∂g₂_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
    #global ∂g₂_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
    #global ∂g₂_∂d = zeros(size(d)) # behövs inte om vi har lokal funktion?
    ∂rᵤ_∂x::Matrix{Float64}
    λᵤ::Vector{Float64}
    λψ::Vector{Float64}
    λᵥₒₗ::Vector{Float64}
end

struct Opt
    m::Int64
    n::Int64
    epsimin::Float64
    xvalue::Vector{Float64}
    xold1::Vector{Float64}
    xold2::Vector{Float64}
    xmin::Vector{Float64}
    xmax::Vector{Float64}
    C::Vector{Float64}
    d2::Vector{Int64}
    a0::Int64
    am::Vector{Int64}
    low::Vector{Float64}
    upp::Vector{Float64}
end



# Initialize the struct with the global variables
