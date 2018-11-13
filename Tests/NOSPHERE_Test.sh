# NOSPHERE_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 04/06/18

# This test compares particle motion in a constant electric and magnetic field without a sphere present in the centre.
# Three cases are investigated: 1) With constant mangetic field 2) constant (central) Electric field
# 3) constant magnetic and electric fields.

# ENSURE THAT TEST_NOSPHERE IS UNCOMMENTED

echo "Test 1:"
echo "Particle Motion conservation in constant magnetic field"
./main -p 0.0 -f 0.0 -m 10 -i 100 -j 100 -se 1 -no 1 > Tests/Data/mag_NOSPHERE_out.txt

echo
echo "Complete!"
echo "Running analytical comparison script Tests/Data/Motion_Test.h"
echo

./Tests/Data/Test -m "BField"

echo
echo "Complete!"
echo "Plotting data from Tests/Data/mag_NOSPHERE_out.txt"
echo

gnuplot -e "filename='Tests/Data/mag_NOSPHERE_out.txt'" -e "prefix='uncharged'" Tests/PlotScripts/NOSPHERE.plt

#echo "Test 2:"
#echo "Particle Motion conservation in constant radial electric field"
#./main -p 0.0 -b 10.0 -m 10 -i 10000 -j 100 -se 1 -no 1 > Tests/Data/elec_NOSPHERE_out.txt
#
#echo
#echo "Complete!"
#echo "Plotting data from Tests/Data/elec_NOSPHERE_out.txt"
#echo
#
#gnuplot -e "filename='Tests/Data/elec_NOSPHERE_out.txt'" -e "prefix='uncharged'" Tests/PlotScripts/NOSPHERE.plt
#
#
#echo
#echo "Test 3:"
#echo "Particle Motion conservation in constant magnetic and electric field"
#
#./main -p -2.5 -b 10.0 -m 10 -u 1.0 -l 1.0 -i 10000 -j 100 -se 1 -no 1 > Tests/Data/magelec_NOSPHERE_out.txt
#
#echo
#echo "Complete!"
#echo "Plotting data from Tests/Data/magelec_NOSPHERE_out.txt"
#echo
#
#gnuplot  -e "filename='Tests/Data/magelec_NOSPHERE_out.txt'" -e "prefix='fixedcharge'" Tests/PlotScripts/NOSPHERE.plt
#
#echo
#echo "Finished!"
#echo 
