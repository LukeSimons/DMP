# ENERGY_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test investigates the degree of energy conservation for the DiMPl code.
# Two cases are investigated: 1) With constant mangetic field 2) constant magnetic and electric fields.
# The distribution of gyro-radii is also plotted

# ENSURE THAT TEST_ENERGY IS UNCOMMENTED

echo "Test 1:"
echo "Energy conservation in constant magnetic field"
./main -p 0.0 -b 10.0 -m 65 -i 10000 -j 100 -se 1 > Tests/Data/mag_ENERGY_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/mag_ENERGY_out.txt"
echo

gnuplot -e "filename='Tests/Data/mag_ENERGY_out.txt'" -e "prefix='uncharged'" Tests/PlotScripts/ENERGY.plt

echo
echo "Test 2:"
echo "Energy conservation in constant magnetic and electric field"

./main -p -2.5 -b 10.0 -m 65 -u 1.0 -l 1.0 -i 10000 -j 100 -se 1 > Tests/Data/magelec_ENERGY_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/magelec_ENERGY_out.txt"
echo

gnuplot  -e "filename='Tests/Data/magelec_ENERGY_out.txt'" -e "prefix='fixedcharge'" Tests/PlotScripts/ENERGY.plt

echo
echo "Test 3:"
echo "Energy conservation in constant magnetic and electric field normalised"

./main -p -2.5 -b 10.0 -m 10 -n 1 -i 10000 -j 100 -se 1 > Tests/Data/magelecnorm_ENERGY_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/magelecnorm_ENERGY_out.txt"
echo

gnuplot  -e "filename='Tests/Data/magelecnorm_ENERGY_out.txt'" -e "prefix='fixedchargenorm'" Tests/PlotScripts/ENERGY.plt

echo
echo "Finished!"
echo 
# -e "set logscale y"
