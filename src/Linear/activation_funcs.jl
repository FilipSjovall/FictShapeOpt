using Plots

# Define the activation functions
sigmoid(x) = 1 / (1 + exp(-x))
relu(x) = max(0, x)
leaky_relu(x, alpha=0.01) = x > 0 ? x : alpha * x
prelu(x, alpha=0.01) = max(alpha * x, x)
elu(x, alpha=1.0) = x > 0 ? x : alpha * (exp(x) - 1)
selu(x, alpha=1.67326, scale=1.0507) = scale * (x > 0 ? x : alpha * (exp(x) - 1))
silu(x) = x * sigmoid(x)
gelu(x) = 0.5x * (1 + tanh(sqrt(2 / π) * (x + 0.044715 * x^3)))
# Generate x values
x = range(-10, 10, length=100)

# Plot the activation functions
plot(
    plot(x, sigmoid.(x), label="Sigmoid",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, tanh.(x), label="Tanh",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, relu.(x), label="ReLU",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, leaky_relu.(x), label="Leaky ReLU",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, prelu.(x), label="PReLU",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, elu.(x), label="ELU",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, selu.(x), label="SELU",xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, silu.(x), label="SiLU", xlims=(-5, 5), ylims=(-1, 1)),
    plot(x, gelu.(x), label="GeLU", xlims=(-5, 5), ylims=(-1, 1)),
    layout=(3, 3), size=(800, 400), framestyle=:origin, legend=:bottomright
)
