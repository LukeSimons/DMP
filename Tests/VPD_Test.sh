# VPD_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test investigates the injected velocity and position distributions for the DiMPl code.
# For ions, electrons and both species together, the injected velocity and position distributions are plotted with gnuplot.
# The distribution of gyro-radii is also plotted

# ENSURE THAT TEST_VELPOSDIST IS UNCOMMENTED IN MAIN.CPP

echo "Test 1:"
echo "Ion velocity distribution"
./main -p 0.0 -b 6.0 -m 100 -i 100000 -j 10000 -c 1.0 > Tests/ions_VPD_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/ions_VPD_out.txt"
echo

gnuplot -e "filename='Tests/ions_VPD_out.txt'" -e "prefix='ion'" Tests/VPD.plt

echo
echo "Test 2:"
echo "Electron velocity distribution"

./main -p 0.0 -b 6.0 -m 100 -i 100000 -j 10000 -c 0.0 > Tests/elecs_VPD_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/elecs_VPD_out.txt"
echo

gnuplot -e "filename='Tests/elecs_VPD_out.txt'" -e "prefix='elec'" -e "bin_width_coeff=50.0" Tests/VPD.plt

echo
echo "Test 3:"
echo "Ion & Electron distribution"
./main -p 0.0 -b 6.0 -m 65 -i 100000 -j 10000 > Tests/VPD_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/VPD_out.txt"
echo

gnuplot -e "set logscale y" -e "filename='Tests/VPD_out.txt'" -e "prefix='both'" -e "bin_width_coeff=50.0" Tests/VPD.plt

echo
echo "Finished!"
echo
