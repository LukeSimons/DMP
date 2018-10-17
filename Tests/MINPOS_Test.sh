# MINPOS_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 11/07/18

# This tests the electric field by comparing the numerical model closest approach value to that of the simulation
# The first test is for ions in constant electric field.

# ENSURE THAT TEST_CLOSEST_APPROACH IS UNCOMMENTED

echo "Test 1:"
echo "Ion minimum positions, constant attractive electric field"

./main -p -0.1 -m 0.0 -z 100.0 -f 5.0 -i 500 -j 500 -c 1.0 -sa 1 -se 1 -no 1 -t 0.0001 > Tests/Data/ions_attract_MINPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/ions_attract_MINPOS_out.txt"
echo

gnuplot -e "filename='Tests/Data/ions_attract_MINPOS_out.txt'" -e "prefix='ion_attract'" Tests/PlotScripts/MINPOS.plt

echo "Test 2:"
echo "Electron minimum positions, constant attractive electric field"

./main -p 0.1 -m 0.0 -z 100.0 -f 5.0 -i 500 -j 500 -c 0.0 -sa 1 -se 1 -no 1 -t 0.0001 > Tests/Data/electrons_attract_MINPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/electrons_attract_MINPOS_out.txt"
echo

gnuplot -e "filename='Tests/Data/electrons_attract_MINPOS_out.txt'" -e "prefix='elec_attract'" Tests/PlotScripts/MINPOS.plt

echo "Test 3:"
echo "Ion minimum positions, constant repulsive electric field"

./main -p 0.1 -m 0.0 -z 100.0 -f 5.0 -i 500 -j 500 -c 1.0 -sa 1 -se 1 -no 1 -t 0.0001 > Tests/Data/ions_repel_MINPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/electrons_repel_MINPOS_out.txt"
echo

gnuplot -e "filename='Tests/Data/electrons_repel_MINPOS_out.txt'" -e "prefix='ion_repel'" Tests/PlotScripts/MINPOS.plt

echo "Test 4:"
echo "Electron minimum positions, constant repulsive electric field"

./main -p -0.1 -m 0.0 -z 100.0 -f 5.0 -i 500 -j 500 -c 0.0 -sa 1 -se 1 -no 1 -t 0.0001 > Tests/Data/electrons_repel_MINPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/electrons_repel_MINPOS_out.txt"
echo

gnuplot -e "filename='Tests/Data/electrons_repel_MINPOS_out.txt'" -e "prefix='elec_repel'" Tests/PlotScripts/MINPOS.plt

echo "Test 5:"
echo "Ion minimum positions, LARGE constant attractive electric field"

./main -p -2.5 -m 0.0 -z 100.0 -f 5.0 -i 500 -j 500 -c 1.0 -sa 1 -se 1 -no 1 -t 0.000001 > Tests/Data/ions_big_attract_MINPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/Data/ions_big_attract_MINPOS_out.txt"
echo

gnuplot -e "filename='Tests/Data/ions_big_attract_MINPOS_out.txt'" -e "prefix='ions_big_attract'" Tests/PlotScripts/MINPOS.plt



echo
echo "Finished!"
echo

