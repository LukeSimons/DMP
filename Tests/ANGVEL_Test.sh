# ENSURE THAT TEST_ANGVEL IS UNCOMMENTED


./main -p 0.0 -b 10.0 -m 65 -i 100000 -j 1000 -se 1 > Tests/ANGVEL_out.txt

gnuplot -e "filename='Tests/ANGVEL_out.txt'" -e "prefix='Avel'" Tests/ANGVEL.plt

