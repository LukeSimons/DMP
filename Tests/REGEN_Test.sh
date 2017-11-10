# ENSURE THAT  TEST_REGEN IS UNCOMMENTED

./main -p 0.0 -b 6.0 -m 100 -i 100000 -j 10000 -c 1.0 > Tests/ions_REGEN_out.txt

gnuplot -e "filename='Tests/ions_REGEN_out.txt'" -e "prefix='ion'" Tests/REGEN.plt

./main -p 0.0 -b 6.0 -m 100 -i 100000 -j 10000 -c 0.0 > Tests/elecs_REGEN_out.txt

gnuplot -e "filename='Tests/elecs_REGEN_out.txt'" -e "prefix='elec'" Tests/REGEN.plt

./main -p 0.0 -b 6.0 -m 65 -i 100000 -j 10000 > Tests/REGEN_out.txt

gnuplot -e "set logscale y" -e "filename='Tests/REGEN_out.txt'" -e "prefix='both'" Tests/REGEN.plt

# High magnetic field test 

./main -p 0.0 -b 6.0 -m 10000 -i 100000 -j 10000 > Tests/REGEN_out.txt

gnuplot -e "filename='Tests/REGEN_out.txt'" -e "prefix='High'" Tests/REGEN.plt


# Low magnetic field test 

./main -p 0.0 -b 12.0 -m 0.1 -i 100000 -j 10000 > Tests/REGEN_out.txt

gnuplot -e "filename='Tests/REGEN_out.txt'" -e "prefix='Low'" Tests/REGEN.plt

