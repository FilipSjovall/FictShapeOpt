using Plots; gr()

x = 0:0.1:π
y1 = x.^2


x_points = [1, 10, 20, 32]
y2 = 7sqrt.(x) 
p=plot(x, y1, color =:steelblue  , label=:none,linewidth=5)
scatter!(x[x_points], y2[x_points], color =:tomato ,label=:none, axis=nothing, marker=:diamond,linewidth=5,markersize=:5)

scatter!(x[x_points], y1[x_points], label=:none, color=:steelblue, marker=:diamond,linewidth=10,markersize=:5)

xl, yl = xlims(p), ylims(p)
annotate!(xl[2],yl[1], text('\u27A4',7,:black,rotation=0))
annotate!(xl[1],yl[2], text('\u27A4',7,:black,rotation=90))
for i in 1:length(x_points)
    plot!([x[x_points[i]], x[x_points[i]]], [y1[x_points[i]], y2[x_points[i]]], color=:red, label=:none,linewidth=2.5, linestyle=:dash)
end
current()  