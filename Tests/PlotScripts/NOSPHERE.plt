if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"

## PLOT DIFFERENCE BETWEEN ANALYTICAL AND NUMERICAL MODEL

set output sprintf("Plots/%s_NoSphere_1.eps",prefix)
# Plot velocity changes as percentage
plot filename using 2
set output sprintf("Plots/%s_NoSphere_2.eps",prefix)
# Plot energy changes as percentage
plot filename using 3

