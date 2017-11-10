if (!exists("filenumi")) filenumi=0
if (!exists("filenumf")) filenumf=100000
if (!exists("iter")) iter=100

filename(n) = sprintf("Data/DMP_Track_%d.txt", n);
#set key off
set logscale x
plot for [i=filenumi:filenumf:iter] filename(i) using 6 with lines lw 1 title sprintf("DMP Track %d.txt",i)
