# Run the program, in this case B=12.22855145T corresponds to Beta=0.3. 
# Running for Ions only

#magfields=("12.18430597" "20.30717662" "28.43004726" "40.61435323" "56.86009452" "81.22870646" "122.8430597" "163.4574129" "203.0717662" "284.3004726" "406.1435323") 
#magfields=("9.747444776" "16.2457413" "22.74403781" "32.49148258" "45.48807562" "64.98296517" "98.27444776" "130.7659303" "162.457413" "227.4403781" "324.9148258")
magfields=("0.3" "0.5" "0.7" "1.0" "1.4" "2.0" "3.0" "4.0" "5.0" "7.0" "10.0")
#magfields=("0.0")
potentials=("0.0")
#potentials=("-0.1" "-0.2" "-0.3" "-0.4" "-0.6" "-0.8" "-1.0" "-1.2" "-1.4" "-1.6" "-1.8" "-2.0" "-2.5" "-3.0" "-3.5" "-4.0" "-5.0" "-7.0" "-9.0" "-11.0" "-13.0" "-15.0" "-17.0" "-19.0" "-22.0" "-25.0" "-30.0" "-45.0" "-60.0" "-80.0" "-100.0" )



for pot in ${potentials[@]}
do
	for BMag in ${magfields[@]}
	do
		./main -p 0.0 -i 250000 -j 250000 -m ${BMag} -n 1 -b 10.0 -o _${BMag}_${pot}.txt -c 1.0 
		cat Data/DMP_${BMag}_${pot}.txt | head -n 37 | tail -n 1 > Tests/Currents.txt
		
#		head -n 37 Data/DMP_${BMag}_${pot}.txt | tail -n 1 > Tests/ParticleCounts.txt
		cat Tests/Currents.txt | awk -v V=${pot} -v bmag=${BMag} '
		{
			print bmag,0.03077729671*bmag,pot, $1,$2;
		}' >> Out.txt
	done
done
