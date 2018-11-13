if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"

## 1 DIMENSIONAL PLOT OF THE ANGULAR VELOCITY IN NORMALISED AND UN-NORMALISED UNITS
set format y "10^{%L}"
set logscale y
set output sprintf("Plots/%s_Norm.eps",prefix)
# Plot Normalised units of Angular Vel
plot filename using 3
set output sprintf("Plots/%s.eps",prefix)
# Plot SI units of Angular vel 
plot filename using 6

