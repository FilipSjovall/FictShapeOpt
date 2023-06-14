fcalls = 0
function f(y, x) # in-place
    global fcalls += 1
    for i in 2:length(x)-1
        y[i] = x[i-1] - 2x[i] + x[i+1]
    end
    y[1] = -2x[1] + x[2]
    y[end] = x[end-1] - 2x[end]
    nothing
end

function g(x) # out-of-place
    global fcalls += 1
    y = zero(x)
    for i in 2:length(x)-1
        y[i] = x[i-1] - 2x[i] + x[i+1]
    end
    y[1] = -2x[1] + x[2]
    y[end] = x[end-1] - 2x[end]
    y
end

using Symbolics
input = rand(30)
output = similar(input)
sparsity_pattern = Symbolics.jacobian_sparsity(f, output, input)
jac = Float64.(sparsity_pattern)

using SparseDiffTools
colors = matrix_colors(jac)
x = rand(30)
forwarddiff_color_jacobian!(jac, f, x, colorvec=colors)
