using LinearAlgebra
N1s(ξ) = (1.0 - ξ) / 2
N2s(ξ) = (1.0 + ξ) / 2
Ns(ξ) = [N1s(ξ) N2s(ξ)]

dN1s(ξ) = -1.0 / 2.0
dN2s(ξ) = 1.0 / 2.0
dNs(ξ) = [dN1s(ξ) dN2s(ξ)]

x1s = [0.9 0.9]
x2s = [1.0 1.0]


xm = [0.95 1.0]

n1 = [-0.7 0.3]
n1 = n1 ./ norm(n1)

n2 = [-0.4 0.6]
n2 = n2 ./ norm(n2)

xs(ξ) = Ns(ξ) * [x1s' x2s']
dxs(ξ) = dNs(ξ) * [x1s' x2s']

ns(ξ) = Ns(ξ) * [n1' n2']
dns(ξ) = dNs(ξ) * [n1' n2']

res(ξ) = cross(vec([xs(ξ) - xm 0.0]), vec([ns(ξ) 0.0]))[3]
dres(ξ) = cross(vec([dxs(ξ) 0.0]), vec([ns(ξ) 0.0]))[3] + cross(vec([xs(ξ) - xm 0.0]), vec([dns(ξ) 0.0]))[3]
res(1.0)
dres(1.0)

ξᵢ = 0.0
ξᵢ₊₁ = 0.0
dξ = 0.0
for i ∈ 1:100
    dξ = -res(ξᵢ) ./ dres(ξᵢ)
    dξ = clamp.(dξ, -0.5, 0.5)

    ξᵢ₊₁ = clamp.(ξᵢ + dξ, -1.0, 1.0)
    println("iteration: $i residual ", res(ξᵢ₊₁), " at point ", ξᵢ₊₁, " and derivative: ", dres(ξᵢ₊₁))
    if norm(res(ξᵢ₊₁)) < 1e-10  #norm(ξᵢ₊₁ - ξᵢ) < 1e-10 && norm(dξ)<1e-6
        println(" norm achieved: ", norm(ξᵢ₊₁ - ξᵢ), " ", ξᵢ₊₁, " ", ξᵢ)
        return ξᵢ₊₁
    end
    ξᵢ = ξᵢ₊₁
end