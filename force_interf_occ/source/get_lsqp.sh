ofn=pep_lsqp.txt

rm -rf ${ofn}

pln=($(grep -ae 'LSQP: '  out_calc.dat | sed -e 's/"//g'| gawk '{print $5}'))
d=($(grep -ae 'LSQP: '  out_calc.dat | sed -e 's/"//g'| gawk '{print $7}'))
echo ${pln}
echo ${d[0]}

len=${#pln[@]}
echo " $len " 

for ((i=0; i<$len ; i++))

do
    vals="$(grep -Eo '[+-]?[0-9]+([.][0-9]+)?' <<<"${pln[i]}")"
    echo ${vals} ${d[i]} >> ${ofn}
done

#charmm gives ax+by+cz=d form for plane equation
# will need to multiply d by 1 in post-processing 
