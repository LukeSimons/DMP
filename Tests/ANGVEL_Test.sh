# ENSURE THAT TEST_ANGVEL IS UNCOMMENTED


./main -p 0.0 -b 6.0 -m 65 -i 1000000 -j 10000 > Tests/ANGVEL_out.txt

gnuplot -e "filename='Tests/ANGVEL_out.txt'" -e "prefix='Avel'" Tests/ANGVEL.plt

