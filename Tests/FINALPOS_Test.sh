# FINALPOS_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test investigates the deviation of the final position to the dust grain surface (which is at unity).
# The first test is for ions in constant magnetic field.
# The second test is for electrons in constant magnetic field.
# The third test is for ions and electrons in constant magnetic field.
# The fourth test is for ions and electrons in a constant magnetic and electric field

# ENSURE THAT TEST_FINALPOS IS UNCOMMENTED

echo "Test 1:"
echo "Ion final positions, constant magnetic field"

./main -p 0.0 -b 10.0 -m 100 -i 10000 -j 1000 -c 1.0 > Tests/ions_FPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/ions_FPOS_out.txt"
echo

gnuplot -e "filename='Tests/ions_FPOS_out.txt'" -e "prefix='ion'" Tests/FINALPOS.plt

echo
echo "Test 2:"
echo "Electron final positions, constant magnetic field"

./main -p 0.0 -b 10.0 -m 100 -i 10000 -j 1000 -c 0.0 > Tests/elecs_FPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/elecs_FPOS_out.txt"
echo

gnuplot -e "filename='Tests/elecs_FPOS_out.txt'" -e "prefix='elec'" Tests/FINALPOS.plt

echo
echo "Test 3:"
echo "Electrons and Ions final positions, constant magnetic field"

./main -p 0.0 -b 10.0 -m 100 -i 10000 -j 1000 > Tests/FPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/FPOS_out.txt"
echo

gnuplot -e "filename='Tests/FPOS_out.txt'" -e "prefix='both'" Tests/FINALPOS.plt

echo
echo "Test 4:"
echo "Electron and Ions final positions, constant magnetic and electric field"

./main -p -2.5 -b 10.0 -m 100 -i 10000 -j 1000 > Tests/charged_FPOS_out.txt

echo
echo "Complete!"
echo "Plotting data from Tests/charged_FPOS_out.txt"
echo

gnuplot -e "filename='Tests/charged_FPOS_out.txt'" -e "prefix='charged_both'" Tests/FINALPOS.plt

echo
echo "Finished!"
echo

