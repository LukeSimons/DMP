# OML_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test investigates the floating potential of dust grains in zero magnetic field.
# The expected normalised potential on the dust grain is ~2.5, or -1750 fundamental electric charges.

# ENSURE THAT TEST_CHARGE &&  IS UNCOMMENTED IN MAIN.CPP

echo "Test 1:"
echo "OML Charging Test"

./main -p 0.0 -i 10000000 -j 1000 -m 0.0 -z 10.0 -f 8.0 > Tests/Data/Out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/Out.txt"
echo

gnuplot Tests/PlotScripts/OML.plt
