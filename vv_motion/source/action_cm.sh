# action_cm.sh: get center of mass of C-domains
# bash action_cm.sh ifn dname

ifn0=${1} 
dname=${2}

fn0=temp_${ifn0}.inp
echo "* ${fn0}" > ${fn0}
echo "*" >> ${fn0}
echo "set rstrm bresi_${ifn0}.str !  assigned from action_cm.sh " >> ${fn0}


if [ "${dname}" = "orig" ]; then # use ../../dcd_nowat
    echo "stream dcd_header.str" >> ${fn0}
else # use ../orie/dcd
    echo "set dname ${dname}" >> ${fn0}
    echo "stream dcd_header_orie.str" >> ${fn0}
fi

echo "stream get_cm.str" >> ${fn0}

charmm  -chsize 100000 < ${fn0} > data_${dname}/out_temp_${ifn0}.dat 

echo "#[xyz]c " > data_${dname}/cm_${ifn0}.dat
grep -a '[XYZ]0 <-\|FRM <-' data_${dname}/out_temp_${ifn0}.dat | \
    sed -e 's/"//g' -e 's/ Parameter: //g' -e 's/<-//g' | \
    gawk '{print $2}' >> data_${dname}/cm_${ifn0}.dat

exit
