if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 1 DIMENSIONAL PLOT OF FINAL POSITIONS

set output sprintf("Plots/%s_FINPOS.eps",prefix)
# Plot ion positions
plot filename using 1

