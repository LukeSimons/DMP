if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 2 DIMENSIONAL PLOTS OF POSITION DISTRIBUTION

set key off

set output sprintf("Plots/%s_EFIELD_Pos_2D.eps",prefix)
# Plot 2D positions
plot filename using 1:2, \
filename using 1:6

set output sprintf("Plots/%s_EFIELD_Vel_2D.eps",prefix)
# Plot 2D positions
plot filename using 1:3, \
filename using 1:9

