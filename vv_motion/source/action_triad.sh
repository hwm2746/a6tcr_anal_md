# action_triad.sh: from coord traj, build triad of V-domains
# bash action_triad.sh dname

dname=${2} # output suffix

#######################################################
ifn1=${1}1
fn0=temp_${ifn1}.inp
echo "* ${fn0}" > ${fn0}
echo "*" >> ${fn0}
echo "set rstrm bresi_${ifn1}.str !  assigned from action_triad.sh " >> ${fn0}

if [ "${dname}" = "orig" ]; then # use ../../dcd_nowat
    echo "stream dcd_header.str" >> ${fn0}
else # use ../orie/dcd
    echo "set dname ${dname}" >> ${fn0}
    echo "stream dcd_header_orie.str" >> ${fn0}
fi

echo "stream get_triad.str" >> ${fn0}

charmm  -chsize 100000 < ${fn0} > data_${dname}/out_temp_${ifn1}.dat

echo "#frame [xyz]d, [xyz]c " > data_${dname}/raw_${ifn1}.dat
grep -a '[XYZ]D <-\|[XYZ]C <-\|FRM <-' data_${dname}/out_temp_${ifn1}.dat  | \
    sed -e 's/"//g' -e 's/ Parameter: //g' -e 's/<-//g' | \
	gawk '{print $2}' >> data_${dname}/raw_${ifn1}.dat

#######################################################
ifn0=${1}0
fn0=temp_${ifn0}.inp
echo "* ${fn0}" > ${fn0}
echo "*" >> ${fn0}
echo "set rstrm bresi_${ifn0}.str !  assigned from action_traid.sh " >> ${fn0}

if [ "${dname}" = "orig" ]; then # use ../../dcd_nowat
    echo "stream dcd_header.str" >> ${fn0}
else # use ../orie/dcd
    echo "set dname ${dname}" >> ${fn0}
    echo "stream dcd_header_orie.str" >> ${fn0}
fi

echo "stream get_triad.str" >> ${fn0}

# don't use '&': grep won't work before finishing charmm
charmm  -chsize 100000 < ${fn0} > data_${dname}/out_temp_${ifn0}.dat

echo "#frame [xyz]d, [xyz]c " > data_${dname}/raw_${ifn0}.dat
grep -a '[XYZ]D <-\|[XYZ]C <-\|FRM <-' data_${dname}/out_temp_${ifn0}.dat  | \
    sed -e 's/"//g' -e 's/ Parameter: //g' -e 's/<-//g' | \
	gawk '{print $2}' >> data_${dname}/raw_${ifn0}.dat
