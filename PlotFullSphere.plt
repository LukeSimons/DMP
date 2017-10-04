if (!exists("filenumi")) filenumi=0
if (!exists("filenumf")) filenumf=5

filename(n) = sprintf("Data/MagPlasData%d.txt", n);
set key off
set parametric
set urange [0:2*pi]
set vrange [0:2*pi]
# Parametric functions for the sphere
r=1
fx(v,u) = r*cos(v)*cos(u)
fy(v,u) = r*cos(v)*sin(u)
fz(v)   = r*sin(v)
splot fx(v,u),fy(v,u),fz(v), \
for [i=filenumi:filenumf] filename(i) using 1:2:3 with lines lw 1

