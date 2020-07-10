#! /bin/bash

# Initial concentration (must agree with the initial condition in OF)
cin=1.0



############################################################################

rm *tmp 2> /dev/null 
rm BT.dat 2> /dev/null

re='[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'

# Produce list of dirs at the various ts

ls | awk -v re=$re '$1 ~ /[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ {print $0}' > ts.tmp

# Sort dirs in numerical orders

sort -g ts.tmp > srt.tmp

# Extract concentration

for n in `cat srt.tmp`
    do
        awk -v n=$n '{if($1 == "internalField"){print n,$3}}'  $n/c >> ff.tmp      
done


# Remove ; at the end of the numbers (OpenFoam feature)

awk -F';' '{ print $1 }' ff.tmp > ff2.tmp

lastc=$(tail -1 ff2.tmp | awk '{print $2}')

# Calculate breaktrough

awk -v ceq=$lastc -v cin=$cin '
    BEGIN{
        s2d=1/3600/24;  # 1/seconds to days 
    }
{ y=(ceq-$2)/(ceq-cin); x=$1*s2d; if(y > 0 && x > 0){print x,y,$2} }' ff2.tmp > BT.dat


rm *tmp 2> /dev/null 
