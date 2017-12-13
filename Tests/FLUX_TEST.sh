# Run the program, in this case B=12.22855145T corresponds to Beta=0.3. 
# Running for Ions only

magfields=("12.22855145" "20.38091908" "28.53328671" "40.76183815" "57.06657341" "81.5236763" "122.2855145" "163.0473526" "203.8091908" "285.3328671" "407.6183815") 
#magfields=("0.0")
potentials=("-0.1" "-0.2" "-0.3" "-0.4" "-0.6" "-0.8" "-1.0" "-1.2" "-1.4" "-1.6" "-1.8" "-2.0" "-2.5" "-3.0" "-3.5" "-4.0" "-5.0" "-7.0" "-9.0" "-11.0" "-13.0" "-15.0" "-17.0" "-19.0" "-22.0" "-25.0" "-30.0" "-45.0" "-60.0" "-80.0" "-100.0" )



for pot in ${potentials[@]}
do
	for BMag in ${magfields[@]}
	do
		./main -p ${pot} -i 1000000 -j 1000000 -m ${BMag} -b 10.0 -o _${BMag}_${pot}.txt -c 1.0
				
		IonIP=$(head -n 20 Data/DMP_${BMag}_${pot}.txt | tail -n 1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')
		IonTemp=$(head -n 18 Data/DMP_${BMag}_${pot}.txt | tail -n 1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')
		
		head -n 37 Data/DMP_${BMag}_${pot}.txt | tail -n 1 > Tests/ParticleCounts.txt
		cat Tests/ParticleCounts.txt | awk -v ip=${IonIP} -v bmag=${BMag} '
		{
			print bmag,0.02453275037*bmag,$1/$9, ip, $1*(ip^2)/(2*$9);
		}' >> Out.txt
	done
done

#			echarge=1.60217662*10^(-19);
#			Mp=1.6726219*10^-27;
#			PI=atan2(0,-1);
#			Flux=2*PI*(ip*10^-6)^(2)*10^(18)*sqrt(echarge/(2*PI*Mp));
#			Current=$1*echarge*Flux/$9;
#			CurrentNorm=4*PI*(10^-6)^(2)*10^18*echarge*sqrt(echarge/(2*PI*Mp));
