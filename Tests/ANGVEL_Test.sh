# FINALPOS_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test Angular velocity of spheres


# ENSURE THAT TEST_ANGVEL IS UNCOMMENTED


./main -p 0.0 -b 10.0 -m 65 -i 100000 -j 1000 -se 1 > Tests/Data/ANGVEL_out.txt

gnuplot -e "filename='Tests/Data/ANGVEL_out.txt'" -e "prefix='Avel'" Tests/PlotScripts/ANGVEL.plt

