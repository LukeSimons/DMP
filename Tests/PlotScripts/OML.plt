set terminal postscript eps size 6.77,5.00 enhanced colour font ",20"

set output "Plots/OML.eps"
plot "Tests/Out.txt" using 1:2
