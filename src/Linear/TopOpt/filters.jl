


function helmholtz(dh,ρ)

    M = assemM(dh)
    K = assemK(dh)
    f = assemFilterLoad(ρ)
    ρ_filtered = solve(K,f)

   return ρ_filtered
end
