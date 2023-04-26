# # # # # # # # # # # # # #
# Run the optimization  # #
# # # # # # # # # # # # # #
g_hist,v_hist,OptIter = Optimize(dh)

# # # # # # # # # 
# Plot history  #
# # # # # # # # # 
using CairoMakie
fig = Figure()
ax1, l1 = lines(fig[1, 1], 1..OptIter, g_hist[1:OptIter], color = :red)
ax2, l2 = lines(fig[2, 1], 1..OptIter, v_hist[1:OptIter], color = :blue)
Legend(fig[1:2, 2], [l1, l2], ["Objective", "Constraint"])
fig