# ENSURE THAT TEST_VELPOSDIST IS UNCOMMENTED

./main -p 0.0 -b 6.0 -m 100 -i 100000 -j 10000 -c 1.0 > Tests/ions_VPD_out.txt

gnuplot -e "filename='Tests/ions_VPD_out.txt'" -e "prefix='ion'" Tests/VPD.plt

./main -p 0.0 -b 6.0 -m 100 -i 100000 -j 10000 -c 0.0 > Tests/elecs_VPD_out.txt

gnuplot -e "filename='Tests/elecs_VPD_out.txt'" -e "prefix='elec'" Tests/VPD.plt

./main -p 0.0 -b 6.0 -m 65 -i 100000 -j 10000 > Tests/VPD_out.txt

gnuplot -e "set logscale y" -e "filename='Tests/VPD_out.txt'" -e "prefix='both'" Tests/VPD.plt
