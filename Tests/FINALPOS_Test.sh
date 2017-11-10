# ENSURE THAT TEST_FINALPOS IS UNCOMMENTED

./main -p 0.0 -b 10.0 -m 100 -i 100000 -j 10000 -c 1.0 > Tests/ions_FPOS_out.txt

gnuplot -e "filename='Tests/ions_FPOS_out.txt'" -e "prefix='ion'" Tests/FINALPOS.plt

./main -p 0.0 -b 10.0 -m 100 -i 100000 -j 10000 -c 0.0 > Tests/elecs_FPOS_out.txt

gnuplot -e "filename='Tests/elecs_FPOS_out.txt'" -e "prefix='elec'" Tests/FINALPOS.plt
