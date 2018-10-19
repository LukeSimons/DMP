if (!exists("filenumi")) filenumi=0
if (!exists("filenumf")) filenumf=5
if (!exists("iter")) iter=1

filename(n) = sprintf("Data/DiMPl_Track_%d.txt", n);
set key off
set parametric
set urange [0:2*pi]
set vrange [0:2*pi]

set xtics font "Helvetica,26"
set ytics font "Helvetica,26"
set ztics font "Helvetica,26"
set xlabel "X" font "Helvetica,30"
set ylabel "Y" font "Helvetica,30"
set zlabel "Z" font "Helvetica,30"

# Parametric functions for the sphere
r=1
fx(v,u) = r*1.0*cos(v)*cos(u)
fy(v,u) = r*1.0*cos(v)*sin(u)
fz(v)   = r*1.0*sin(v)
splot fx(v,u),fy(v,u),fz(v), \
for [i=filenumi:filenumf:iter] filename(i) using 1:2:3 
#for [i=filenumi:filenumf:iter] filename(i) using 1:2:3 with lines lw 2




#set xrange [-11:11]
#set yrange [-11:11]
#set zrange [-11:11]
replot
