if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 1 DIMENSIONAL PLOT OF THERMAL CURRENT

# Each bar is half the (visual) width of its x-range.
set boxwidth 0.005 absolute
set style fill solid 1.0 noborder
bin_width = 0.005;
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

set output sprintf("Plots/%s_FINPOS.eps",prefix)
# Plot ion positions
plot filename using (rounded($1)):(1) smooth frequency with boxes

