if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"
if (!exist("bin_width_coeff")) bin_width_coeff=1.0
set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF POSITION AND VELOCITY DISTRIBUTIONS

set output sprintf("Plots/%s_VPD_Pos.eps",prefix)
# Plot positions
splot filename using 1:2:3

set output sprintf("Plots/%s_VPD_Vel.eps",prefix)
# Plot Velocities
splot filename using 4:5:6


## 2 DIMENSIONAL PLOTS OF POSITION DISTRIBUTION

set output sprintf("Plots/%s_VPD_Pos_2D.eps",prefix)
# Plot 2D positions
plot filename using 1:2


## 1 DIMENSIONAL PLOTS OF POSITION AND VELOCITY DISTRIBUTION

# Each bar is half the (visual) width of its x-range.
set boxwidth 0.05 absolute
set style fill solid 1.0 noborder
bin_width=0.1
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

set output sprintf("Plots/%s_VPD_Pos_x.eps",prefix)
# Plot histogram of x positions
plot filename using (rounded($1)):(1) smooth frequency with boxes
set output sprintf("Plots/%s_VPD_Pos_y.eps",prefix)
# Plot histogram of y positions
plot filename using (rounded($2)):(2) smooth frequency with boxes
set output sprintf("Plots/%s_VPD_Pos_z.eps",prefix)
# Plot histogram of z positions
plot filename using (rounded($3)):(3) smooth frequency with boxes

set boxwidth 0.05*bin_width_coeff absolute
bin_width=0.1*bin_width_coeff
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

set output sprintf("Plots/%s_VPD_Vel_x.eps",prefix)
# Plot histogram of x Velocities
plot filename using (rounded($4)):(4) smooth frequency with boxes
set output sprintf("Plots/%s_VPD_Vel_y.eps",prefix)
# Plot histogram of y Velocities
plot filename using (rounded($5)):(5) smooth frequency with boxes
set output sprintf("Plots/%s_VPD_Vel_z.eps",prefix)
# Plot histogram of z Velocities
plot filename using (rounded($6)):(6) smooth frequency with boxes


set boxwidth 0.05/bin_width_coeff absolute
bin_width=0.1/bin_width_coeff
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

set output sprintf("Plots/%s_VPD_GyroRad.eps",prefix)
# Plot histogram of Gyro-radii
plot filename using (rounded($7)):(7) smooth frequency with boxes
