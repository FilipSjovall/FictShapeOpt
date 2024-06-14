function simp(ρ,q)
    return ρ^q
end

function d_simp(ρ,q)
    return q*ρ^(q-1)
end

function compute_moduli(E,ν)
    G = E/(2(1+ν))
    K = E/(3(1-2ν))
    return K,G
end

function ramp(x,q,δ₀)
    Χ = δ₀ + (1-δ₀)*x/(1+q*(1-x))
    return Χ
end

function d_ramp(x,q,δ₀)
    # ...
end

function heaviside(x,β,η)
    H = ( tanh(β*η) + tanh(β*(x-η)) ) / ( tanh(β*η) + tanh(β*(1-η)) )
    return H
end

function d_heaviside(x,β,η)
    # ...
end
