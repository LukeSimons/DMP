if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"

## 1 DIMENSIONAL PLOT OF PERCENTAGE UNCERTAINTY IN ENERGY

set output sprintf("Plots/%s_Energy_Kin.eps",prefix)
# Plot velocity changes as percentage
plot filename using 2
set output sprintf("Plots/%s_Energy_Tot.eps",prefix)
# Plot energy changes as percentage
plot filename using 3

