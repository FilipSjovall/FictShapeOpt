### Uppgift 7.3
if 1 == 1
    ω = 100 * 2 * π
    D = 200e-3
    d = 100e-3
    t = 20e-3
    ν = 0.29
    ρ = 7850

    Rot(x) =-(3+ν)/8 * ρ * ω^2 * x^2/4
    Rot2(x)=-(1+3ν)/8 * ρ * ω^2 * x^2/4
    #Rot(d)
    #Rot(D)

    B = ( Rot(d) - Rot(D) ) / (4/d^2 - 4/D^2)

    A = 4B/d^2 - Rot(d)
    A = 4B/D^2 - Rot(D)

    σᵩ(x) = Rot2(x) + A + 4B/x^2


end

σᵩ(d)
σᵩ(D)
