# ENERGY_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test investigates the degree of energy conservation for the DiMPl code.
# Two cases are investigated: 1) With constant mangetic field 2) constant magnetic and electric fields.
# The distribution of gyro-radii is also plotted

# ENSURE THAT TEST_ENERGY IS UNCOMMENTED

echo "Test 1:"
echo "Energy conservation in constant magnetic field"
./main -p 0.0 -b 10.0 -m 65 -i 10000 -j 100 -se 1 > Tests/mag_ENERGY_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/mag_ENERGY_out.txt"
echo

gnuplot -e "filename='Tests/mag_ENERGY_out.txt'" -e "prefix='uncharged'" Tests/ENERGY.plt

echo
echo "Test 2:"
echo "Energy conservation in constant magnetic and electric field"

./main -p -2.5 -b 10.0 -m 65 -u 1.0 -l 1.0 -i 10000 -j 100 -se 1 > Tests/magelec_ENERGY_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/magelec_ENERGY_out.txt"
echo

gnuplot  -e "filename='Tests/magelec_ENERGY_out.txt'" -e "prefix='fixedcharge'" Tests/ENERGY.plt

echo
echo "Finished!"
echo 
# -e "set logscale y"
