if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF POSITION AND VELOCITY DISTRIBUTIONS

set output sprintf("Plots/%s_VPD_Pos.eps",prefix)
# Plot ion positions
splot filename using 1:2:3

set output sprintf("Plots/%s_VPD_Vel.eps",prefix)
# Plot ion Velocities
splot filename using 4:5:6


## 2 DIMENSIONAL PLOTS OF POSITION DISTRIBUTION

set output sprintf("Plots/%s_VPD_Pos_2D.eps",prefix)
# Plot ion positions
plot filename using 1:2


## 1 DIMENSIONAL PLOTS OF POSITION AND VELOCITY DISTRIBUTION

set output sprintf("Plots/%s_VPD_Pos_x.eps",prefix)
# Plot ion positions
plot filename using 1
set output sprintf("Plots/%s_VPD_Pos_y.eps",prefix)
# Plot ion positions
plot filename using 2
set output sprintf("Plots/%s_VPD_Pos_z.eps",prefix)
# Plot ion positions
plot filename using 3

set output sprintf("Plots/%s_VPD_Vel_x.eps",prefix)
# Plot ion Velocities
plot filename using 4
set output sprintf("Plots/%s_VPD_Vel_y.eps",prefix)
# Plot ion Velocities
plot filename using 5
set output sprintf("Plots/%s_VPD_Vel_z.eps",prefix)
# Plot ion Velocities
plot filename using 6

set output sprintf("Plots/%s_VPD_GyroRad.eps",prefix)
# Plot Gyro-radii
plot filename using 7
