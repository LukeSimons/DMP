# ENSURE THAT TEST_ENERGY IS UNCOMMENTED


./main -p 0.0 -b 6.0 -m 65 -i 1000000 -j 10000 > Tests/ENERGY_out.txt

gnuplot -e "filename='Tests/ENERGY_out.txt'" -e "prefix='uncharged'" Tests/ENERGY.plt

./main -p -2.5 -b 6.0 -m 65 -i 1000000 -j 10000 > Tests/ENERGY_out.txt
gnuplot  -e "filename='Tests/ENERGY_out.txt'" -e "prefix='fixedcharge'" Tests/ENERGY.plt
# -e "set logscale y"
