# EFIELD_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 18/06/18

# This test compares the motion of particles in DiMPle in an Electric field with their expected analytical motion.
# This is used as a validation test that the code is behaving as expected

# ENSURE THAT POINT_INJECTION && SAVE_TRACKS ARE UNCOMMENTED

echo
echo "Test motion of initially stationary particle in electric field, compare to Euler Numerical integration"
echo 

./main -p -1.0 -m 0.0 -z 5.0 -i 1 -j 1 -c 1.0 -sa 1

echo
echo "Complete!"
echo

IonTimestepLine=$(head -n 24 Data/DiMPl.txt | tail -n 1)

IonTimeStep=$(echo ${IonTimestepLine:13:10} | sed -e 's/[eE]+*/*10^/' | bc -l )
ChargesOnSphere=$(head -n 31 Data/DiMPl.txt | tail -n 1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')


./test -m "EField" -c $ChargesOnSphere -t $IonTimeStep > Tests/Data/ions_EFIELD.txt


gnuplot -e "filename='Tests/Data/ions_EFIELD.txt'" -e "prefix='NoDrift'" Tests/PlotScripts/EFIELD.plt

echo
echo "Finished!"
echo


echo
echo "Test motion of moving particle in electric field, compare to Euler Numerical integration"
echo 

./main -p -1.0 -m 0.0 -z 5.0 -v -1000.0 -i 1 -j 1 -c 1.0 -sa 1

echo
echo "Complete!"
echo

IonTimestepLine=$(head -n 24 Data/DiMPl.txt | tail -n 1)

IonTimeStep=$(echo ${IonTimestepLine:13:10} | sed -e 's/[eE]+*/*10^/' | bc -l )
ChargesOnSphere=$(head -n 31 Data/DiMPl.txt | tail -n 1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')


./test -m "EField" -c $ChargesOnSphere -t $IonTimeStep > Tests/Data/ions_EFIELD_drift_out.txt


gnuplot -e "filename='Tests/Data/ions_EFIELD_drift_out.txt'" -e "prefix='Drift'" Tests/PlotScripts/EFIELD.plt

echo
echo "Finished!"
echo

