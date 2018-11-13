# THERMALCURRENT_Test.sh
# Author: Luke Simons (ls5115@ic.ac.uk)
# Date: 25/04/18

# This test Calculates the current to an uncharged sphere outside of a magnetic field.
# This is used as a validation test that the code is behaving as expected

echo "Test 1:"
echo "Thermal Current to uncharged sphere, close to sphere"
echo 

echo "Warning! This test will take an extremely long time to complete!"

./main -p 0.0 -m 0.0 -z 0.5 -f 5.0 -i 5000000 -j 40000 -c 1.0 -sa 1 -se 1

echo
echo "Complete!"
echo "Plotting data from Data/DiMPl.txt"
echo

head -n 48 Data/DiMPl.txt | tail -n 1 > Tests/Data/CurrentData.txt
head -n 55 Data/DiMPl.txt | tail -n 1 > Tests/Data/CountingData.txt

Current="$(paste Tests/Data/CurrentData.txt | awk '
{
print $1;
}' )"
CollectedNumErr="$(paste Tests/Data/CountingData.txt | awk '
{
print 1/sqrt($1);
}' )"
paste Tests/Data/CountingData.txt | awk '
{
print 1/sqrt($1);
}'> Out.txt
echo
echo "Thermal Current = $Current +/- $CollectedNumErr"
echo

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

UpperLim=1.01
LowerLim=0.99

Bool1=$(echo $Current'>'$LowerLim | bc -l)
Bool2=$(echo $Current'<'$UpperLim | bc -l)

echo
echo "Finished!"
echo

if (( $(echo $Current'<'$LowerLim | bc -l) ));
then
	echo -e "${RED}Test Failed!${NC}"
elif(( $(echo $Current'>'$UpperLim | bc -l) ));
then
	echo -e "${RED}Test Failed!${NC}"
else
	echo -e "${GREEN}Test Sucessful!${NC}"
fi;

echo 
echo "Second Test"
echo "Thermal Current ot uncharged sphere, far from sphere"
echo

./main -p 0.0 -m 0.0 -z 5.0 -f 10.0 -i 5000000 -j 40000 -c 1.0 -sa 1 -se 1

echo
echo "Complete!"
echo "Plotting data from Data/DiMPl.txt"
echo

head -n 48 Data/DiMPl.txt | tail -n 1 > Tests/Data/CurrentData.txt
head -n 55 Data/DiMPl.txt | tail -n 1 > Tests/Data/CountingData.txt

Current="$(paste Tests/Data/CurrentData.txt | awk '
{
print $1;
}' )"
CollectedNumErr="$(paste Tests/Data/CountingData.txt | awk '
{
print 1/sqrt($1);
}' )"
paste Tests/Data/CountingData.txt | awk '
{
print 1/sqrt($1);
}'> Out.txt
echo
echo "Thermal Current = $Current +/- $CollectedNumErr"
echo

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

UpperLim=1.01
LowerLim=0.99

Bool1=$(echo $Current'>'$LowerLim | bc -l)
Bool2=$(echo $Current'<'$UpperLim | bc -l)

echo
echo "Finished!"
echo

if (( $(echo $Current'<'$LowerLim | bc -l) ));
then
	echo -e "${RED}Test Failed!${NC}"
elif(( $(echo $Current'>'$UpperLim | bc -l) ));
then
	echo -e "${RED}Test Failed!${NC}"
else
	echo -e "${GREEN}Test Sucessful!${NC}"
fi;

