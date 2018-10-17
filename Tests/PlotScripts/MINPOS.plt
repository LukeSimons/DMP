if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"

set output sprintf("Plots/%s_MINPOS.eps",prefix)
## 1 DIMENSIONAL PLOT OF MINIMUM POSITIONS

set ylabel("Normalised Deviation")
set xlabel("Index of Particle")
plot filename using ($3-($4*1e6))/$3

