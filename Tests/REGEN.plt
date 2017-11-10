if (!exist("filename")) filename="NAN.txt"
if (!exist("prefix")) prefix="NAN"

set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"


## 3 DIMENSIONAL PLOTS OF POSITION AND VELOCITY DISTRIBUTIONS

set output sprintf("Plots/%s_Regen_Pos_Regen.eps",prefix)
# Plot post regeneration ion positions
splot filename using 1:2:3

set output sprintf("Plots/%s_Regen_Vel_Regen.eps",prefix)
# Plot post regeneration ion Velocities
splot filename using 4:5:6


## 2 DIMENSIONAL PLOTS OF POSITION DISTRIBUTION

set output sprintf("Plots/%s_Regen_Pos_2D.eps",prefix)
# Plot post regeneration ion positions
plot filename using 1:2


## 1 DIMENSIONAL PLOTS OF POSITION AND VELOCITY DISTRIBUTION

set output sprintf("Plots/%s_Regen_Pos_x.eps",prefix)
# Plot post regeneration ion positions
plot filename using 1
set output sprintf("Plots/%s_Regen_Pos_y.eps",prefix)
# Plot post regeneration ion positions
plot filename using 2
set output sprintf("Plots/%s_Regen_Pos_z.eps",prefix)
# Plot post regeneration ion positions
plot filename using 3

set output sprintf("Plots/%s_Regen_Vel_x.eps",prefix)
# Plot post regeneration ion Velocities
plot filename using 4
set output sprintf("Plots/%s_Regen_Vel_y.eps",prefix)
# Plot post regeneration ion Velocities
plot filename using 5
set output sprintf("Plots/%s_Regen_Vel_z.eps",prefix)
# Plot post regeneration ion Velocities
plot filename using 6

set output sprintf("Plots/%s_Regen_GyroRad.eps",prefix)
# Plot post regeneration gyro-radii
plot filename using 7
