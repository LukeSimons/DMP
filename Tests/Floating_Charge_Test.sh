# Insure that TEST_CHARGING and SELF_CONS_CHARGE have been defined

#./main -p 0.0 -b 10.0 -m 50 -i 100000 -j 10000 -d 0.000001 > Charge_out.txt

gnuplot Tests/PlotScripts/PlotZvelocities.plt
gnuplot plot "Charge_out.txt" using 1:2
